#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <windowsx.h>
#include "motor_synergy.h"
#include <vector>
#include <string>
#include <fstream>
#include <format>
#include <cmath>

constexpr COLORREF COL_BG       = RGB(22, 22, 28);
constexpr COLORREF COL_GRID     = RGB(40, 40, 50);
constexpr COLORREF COL_ALGO     = RGB(100, 149, 237);
constexpr COLORREF COL_WIND     = RGB(255, 99, 71);
constexpr COLORREF COL_HUMAN    = RGB(34, 200, 80);
constexpr COLORREF COL_START    = RGB(0, 255, 100);
constexpr COLORREF COL_TARGET   = RGB(255, 165, 0);
constexpr COLORREF COL_TEXT     = RGB(190, 190, 200);
constexpr COLORREF COL_TEXTDIM  = RGB(120, 120, 135);
constexpr COLORREF COL_GRAPH_BG = RGB(18, 18, 22);

constexpr int TIMER_ANIM = 1;
constexpr int GRAPH_H = 160;
constexpr int STATUS_H = 52;

struct state {
    double sx = 150, sy = 280;
    double ex = 780, ey = 280;
    double tw = 20.0;

    std::vector<motor_synergy::trajectory_point> algo, wind, human;
    motor_synergy::metrics algo_m{}, wind_m{};

    bool recording = false;
    LARGE_INTEGER rec_start{}, perf_freq{};

    bool animating = false;
    size_t anim_idx = 0;
    DWORD anim_start = 0;

    motor_synergy::config cfg;
};
static state g;
static HWND g_hwnd;
static HFONT g_font, g_font_sm;

struct speed_sample { double t, v; };

std::vector<speed_sample> compute_speeds(const std::vector<motor_synergy::trajectory_point>& p) {
    std::vector<speed_sample> out;
    for (size_t i = 1; i < p.size(); ++i) {
        double dx = p[i].x - p[i - 1].x, dy = p[i].y - p[i - 1].y;
        double dt = p[i].t - p[i - 1].t;
        if (dt > 0.0) out.push_back({p[i].t, std::hypot(dx, dy) / dt});
    }
    return out;
}

void fill_rect(HDC dc, int x, int y, int w, int h, COLORREF c) {
    HBRUSH br = CreateSolidBrush(c);
    RECT r = {x, y, x + w, y + h};
    FillRect(dc, &r, br);
    DeleteObject(br);
}

void draw_circle(HDC dc, int cx, int cy, int r, COLORREF col, bool filled) {
    HPEN pen = CreatePen(PS_SOLID, 2, col);
    HBRUSH br = filled ? CreateSolidBrush(col) : (HBRUSH)GetStockObject(NULL_BRUSH);
    SelectObject(dc, pen);
    SelectObject(dc, br);
    Ellipse(dc, cx - r, cy - r, cx + r, cy + r);
    DeleteObject(pen);
    if (filled) DeleteObject(br);
}

void draw_path(HDC dc, const std::vector<motor_synergy::trajectory_point>& pts,
               COLORREF col, size_t count, int canvas_h)
{
    if (pts.empty() || count == 0) return;
    size_t n = std::min(count, pts.size());

    HPEN pen = CreatePen(PS_SOLID, 2, col);
    HPEN old_pen = (HPEN)SelectObject(dc, pen);

    MoveToEx(dc, (int)pts[0].x, (int)pts[0].y, nullptr);
    for (size_t i = 1; i < n; ++i)
        LineTo(dc, (int)pts[i].x, (int)pts[i].y);

    HBRUSH dot = CreateSolidBrush(col);
    SelectObject(dc, dot);
    for (size_t i = 0; i < n; ++i) {
        int px = (int)pts[i].x, py = (int)pts[i].y;
        Ellipse(dc, px - 2, py - 2, px + 2, py + 2);
    }

    if (n > 0 && n < pts.size()) {
        int hx = (int)pts[n - 1].x, hy = (int)pts[n - 1].y;
        HBRUSH head = CreateSolidBrush(RGB(255, 255, 255));
        SelectObject(dc, head);
        Ellipse(dc, hx - 5, hy - 5, hx + 5, hy + 5);
        DeleteObject(head);
    }

    SelectObject(dc, old_pen);
    DeleteObject(pen);
    DeleteObject(dot);
}

void draw_speed_graph(HDC dc, int gx, int gy, int gw, int gh,
    const std::vector<speed_sample>& s_algo,
    const std::vector<speed_sample>& s_wind,
    const std::vector<speed_sample>& s_human)
{
    fill_rect(dc, gx, gy, gw, gh, COL_GRAPH_BG);

    double max_t = 1.0, max_v = 0.01;
    auto update_max = [&](const std::vector<speed_sample>& s) {
        for (auto& p : s) {
            if (p.t > max_t) max_t = p.t;
            if (p.v > max_v) max_v = p.v;
        }
    };
    update_max(s_algo);
    update_max(s_wind);
    update_max(s_human);

    SelectObject(dc, g_font_sm);
    SetTextColor(dc, COL_TEXTDIM);
    SetBkMode(dc, TRANSPARENT);

    auto label = std::format(L"peak: {:.2f} px/ms", max_v);
    TextOutW(dc, gx + 4, gy + 2, label.c_str(), (int)label.size());

    auto label_t = std::format(L"time: {:.0f} ms", max_t);
    TextOutW(dc, gx + gw - 120, gy + gh - 16, label_t.c_str(), (int)label_t.size());

    HPEN grid_pen = CreatePen(PS_DOT, 1, COL_GRID);
    SelectObject(dc, grid_pen);
    for (int i = 1; i <= 3; ++i) {
        int y = gy + gh - (int)(gh * i / 4.0);
        MoveToEx(dc, gx, y, nullptr);
        LineTo(dc, gx + gw, y);
    }
    DeleteObject(grid_pen);

    auto draw_curve = [&](const std::vector<speed_sample>& samples, COLORREF col) {
        if (samples.size() < 2) return;
        HPEN pen = CreatePen(PS_SOLID, 2, col);
        SelectObject(dc, pen);
        bool first = true;
        for (auto& s : samples) {
            int px = gx + (int)(s.t / max_t * gw);
            int py = gy + gh - (int)(s.v / max_v * (gh - 8));
            py = std::clamp(py, gy, gy + gh);
            if (first) { MoveToEx(dc, px, py, nullptr); first = false; }
            else LineTo(dc, px, py);
        }
        DeleteObject(pen);
    };

    draw_curve(s_algo, COL_ALGO);
    draw_curve(s_wind, COL_WIND);
    draw_curve(s_human, COL_HUMAN);
}

void text_line(HDC dc, int x, int& y, const std::wstring& s, COLORREF col = COL_TEXT) {
    SetTextColor(dc, col);
    TextOutW(dc, x, y, s.c_str(), (int)s.size());
    y += 16;
}

void draw_all(HDC dc, int w, int h) {
    fill_rect(dc, 0, 0, w, h, COL_BG);

    int canvas_h = h - GRAPH_H - STATUS_H;

    HPEN grid_pen = CreatePen(PS_DOT, 1, COL_GRID);
    SelectObject(dc, grid_pen);
    for (int x = 0; x < w; x += 100) { MoveToEx(dc, x, 0, nullptr); LineTo(dc, x, canvas_h); }
    for (int y = 0; y < canvas_h; y += 100) { MoveToEx(dc, 0, y, nullptr); LineTo(dc, w, y); }
    DeleteObject(grid_pen);

    HPEN div_pen = CreatePen(PS_SOLID, 1, RGB(60, 60, 70));
    SelectObject(dc, div_pen);
    MoveToEx(dc, 0, canvas_h, nullptr); LineTo(dc, w, canvas_h);
    DeleteObject(div_pen);

    draw_circle(dc, (int)g.sx, (int)g.sy, 6, COL_START, true);
    draw_circle(dc, (int)g.ex, (int)g.ey, (int)g.tw, COL_TARGET, false);

    size_t algo_count = g.animating ? g.anim_idx : g.algo.size();
    draw_path(dc, g.algo, COL_ALGO, algo_count, canvas_h);
    draw_path(dc, g.wind, COL_WIND, g.wind.size(), canvas_h);
    draw_path(dc, g.human, COL_HUMAN, g.human.size(), canvas_h);

    auto sp_a = compute_speeds(g.algo);
    auto sp_w = compute_speeds(g.wind);
    auto sp_h = compute_speeds(g.human);
    draw_speed_graph(dc, 0, canvas_h + 1, w, GRAPH_H, sp_a, sp_w, sp_h);

    int ty = canvas_h + GRAPH_H + 4;
    SelectObject(dc, g_font_sm);
    SetBkMode(dc, TRANSPARENT);

    if (!g.algo.empty()) {
        auto& m = g.algo_m;
        auto s = std::format(L"[SigmaDrift]  MT={:.0f}ms  Fitts={:.0f}ms  PL={:.3f}  PkSpd={:.2f}px/ms  SubMov={}  EndErr={:.1f}px",
            m.movement_time, m.fitts_predicted_mt, m.path_efficiency, m.peak_speed, m.num_submovements, m.endpoint_error);
        text_line(dc, 8, ty, s, COL_ALGO);
    }
    if (!g.wind.empty()) {
        auto& m = g.wind_m;
        auto s = std::format(L"[WindMouse]   MT={:.0f}ms  PL={:.3f}  PkSpd={:.2f}px/ms  SubMov={}  EndErr={:.1f}px",
            m.movement_time, m.path_efficiency, m.peak_speed, m.num_submovements, m.endpoint_error);
        text_line(dc, 8, ty, s, COL_WIND);
    }

    auto help = std::format(L"LClick=Start  RClick=Target  Space=Generate  W=WindMouse  R=Record  C=Clear  S=CSV  +/-=Width({:.0f})", g.tw);
    text_line(dc, 8, ty, help, COL_TEXTDIM);

    if (g.recording) {
        SetTextColor(dc, RGB(255, 60, 60));
        TextOutW(dc, w - 140, canvas_h + GRAPH_H + 4, L"* RECORDING *", 13);
    }
}

void save_csv() {
    auto write = [](const char* name, const std::vector<motor_synergy::trajectory_point>& pts) {
        std::string path = std::string("E:\\Code\\MouseMovementAlgo\\") + name; // change this to your path if you want exports
        std::ofstream f(path);
        f << "t_ms,x,y\n";
        for (auto& p : pts)
            f << std::format("{:.3f},{:.4f},{:.4f}\n", p.t, p.x, p.y);
    };
    if (!g.algo.empty()) write("trajectory_sigmadrift.csv", g.algo);
    if (!g.wind.empty()) write("trajectory_windmouse.csv", g.wind);
    if (!g.human.empty()) write("trajectory_human.csv", g.human);
}

LRESULT CALLBACK wndproc(HWND hw, UINT msg, WPARAM wp, LPARAM lp) {
    switch (msg) {
    case WM_ERASEBKGND: return 1;

    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hw, &ps);
        RECT rc; GetClientRect(hw, &rc);
        int w = rc.right, h = rc.bottom;
        HDC mem = CreateCompatibleDC(hdc);
        HBITMAP bmp = CreateCompatibleBitmap(hdc, w, h);
        HBITMAP old = (HBITMAP)SelectObject(mem, bmp);
        draw_all(mem, w, h);
        BitBlt(hdc, 0, 0, w, h, mem, 0, 0, SRCCOPY);
        SelectObject(mem, old);
        DeleteObject(bmp);
        DeleteDC(mem);
        EndPaint(hw, &ps);
        return 0;
    }

    case WM_LBUTTONDOWN:
        if (!g.recording) {
            g.sx = GET_X_LPARAM(lp); g.sy = GET_Y_LPARAM(lp);
            InvalidateRect(hw, nullptr, FALSE);
        }
        return 0;

    case WM_RBUTTONDOWN:
        g.ex = GET_X_LPARAM(lp); g.ey = GET_Y_LPARAM(lp);
        InvalidateRect(hw, nullptr, FALSE);
        return 0;

    case WM_MOUSEMOVE:
        if (g.recording) {
            LARGE_INTEGER now;
            QueryPerformanceCounter(&now);
            double ms = (double)(now.QuadPart - g.rec_start.QuadPart) / g.perf_freq.QuadPart * 1000.0;
            g.human.push_back({(double)GET_X_LPARAM(lp), (double)GET_Y_LPARAM(lp), ms});
            InvalidateRect(hw, nullptr, FALSE);
        }
        return 0;

    case WM_TIMER:
        if (wp == TIMER_ANIM && g.animating) {
            DWORD elapsed = GetTickCount() - g.anim_start;
            while (g.anim_idx < g.algo.size() && g.algo[g.anim_idx].t <= (double)elapsed)
                ++g.anim_idx;
            if (g.anim_idx >= g.algo.size()) {
                g.animating = false;
                KillTimer(hw, TIMER_ANIM);
            }
            InvalidateRect(hw, nullptr, FALSE);
        }
        return 0;

    case WM_KEYDOWN:
        switch (wp) {
        case VK_SPACE: {
            double dist = std::hypot(g.ex - g.sx, g.ey - g.sy);
            g.cfg.target_width = g.tw;
            g.algo = motor_synergy::generate(g.sx, g.sy, g.ex, g.ey, g.cfg);
            g.algo_m = motor_synergy::compute_metrics(g.algo, g.ex, g.ey, g.tw, dist);
            g.anim_idx = 0;
            g.animating = true;
            g.anim_start = GetTickCount();
            SetTimer(hw, TIMER_ANIM, 16, nullptr);
            InvalidateRect(hw, nullptr, FALSE);
            break;
        }
        case 'W': {
            double dist = std::hypot(g.ex - g.sx, g.ey - g.sy);
            g.wind = windmouse::generate(g.sx, g.sy, g.ex, g.ey);
            g.wind_m = motor_synergy::compute_metrics(g.wind, g.ex, g.ey, g.tw, dist);
            InvalidateRect(hw, nullptr, FALSE);
            break;
        }
        case 'R':
            g.recording = !g.recording;
            if (g.recording) {
                g.human.clear();
                QueryPerformanceFrequency(&g.perf_freq);
                QueryPerformanceCounter(&g.rec_start);
            }
            InvalidateRect(hw, nullptr, FALSE);
            break;
        case 'C':
            g.algo.clear(); g.wind.clear(); g.human.clear();
            g.animating = false;
            KillTimer(hw, TIMER_ANIM);
            InvalidateRect(hw, nullptr, FALSE);
            break;
        case 'S':
            save_csv();
            MessageBoxW(hw, L"Trajectories exported to CSV.", L"Saved", MB_OK);
            break;
        case VK_OEM_PLUS: case VK_ADD:
            g.tw = std::min(g.tw + 2.0, 100.0);
            InvalidateRect(hw, nullptr, FALSE);
            break;
        case VK_OEM_MINUS: case VK_SUBTRACT:
            g.tw = std::max(g.tw - 2.0, 4.0);
            InvalidateRect(hw, nullptr, FALSE);
            break;
        case VK_ESCAPE:
            PostQuitMessage(0);
            break;
        }
        return 0;

    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcW(hw, msg, wp, lp);
}

int WINAPI wWinMain(HINSTANCE inst, HINSTANCE, LPWSTR, int show) {
    g_font    = CreateFontW(16, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, 0, 0, CLEARTYPE_QUALITY, FIXED_PITCH, L"Consolas");
    g_font_sm = CreateFontW(14, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, 0, 0, CLEARTYPE_QUALITY, FIXED_PITCH, L"Consolas");

    WNDCLASSEXW wc{};
    wc.cbSize = sizeof(wc);
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = wndproc;
    wc.hInstance = inst;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.lpszClassName = L"SigmaDriftHarness";
    RegisterClassExW(&wc);

    g_hwnd = CreateWindowExW(0, wc.lpszClassName, L"SigmaDrift - Harness",
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 1050, 720,
        nullptr, nullptr, inst, nullptr);

    ShowWindow(g_hwnd, show);
    UpdateWindow(g_hwnd);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }

    DeleteObject(g_font);
    DeleteObject(g_font_sm);
    return (int)msg.wParam;
}

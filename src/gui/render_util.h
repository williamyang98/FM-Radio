#pragma once

#include <complex>
#include <imgui.h>
#include <implot.h>
#include <fmt/core.h>
#include "utility/span.h"

static void PlotComplexLine(const char *label, const std::complex<float>* x, const int N) {
    auto* x0 = reinterpret_cast<const float*>(x);
    auto I_label = fmt::format("{:s} (I)", label);
    auto Q_label = fmt::format("{:s} (Q)", label);
    ImPlot::PlotLine(I_label.c_str(), &x0[0], N, 1, 0, 0, 0, sizeof(std::complex<float>));
    ImPlot::PlotLine(Q_label.c_str(), &x0[1], N, 1, 0, 0, 0, sizeof(std::complex<float>));
}

static void PlotConstellationSubplot(tcb::span<const std::complex<float>> x, const char* label, const double A=4) {
    const int TOTAL_ROWS = 1;
    const int TOTAL_COLS = 2;
    static double y_min = -A;
    static double y_max = A;

    auto* buf = reinterpret_cast<const float*>(x.data());
    const int N = (int)x.size();

    if (ImPlot::BeginSubplots(label, TOTAL_ROWS, TOTAL_COLS, ImVec2(-1,-1))) {
        if (ImPlot::BeginPlot("Constellation", ImVec2(-1,0), ImPlotFlags_Equal)) {
            ImPlot::SetupAxisLimits(ImAxis_X1, -A, A, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -A, A, ImPlotCond_Once);
            ImPlot::SetupAxisLinks(ImAxis_Y1, &y_min, &y_max);
            const float marker_size = 2.0f;
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Cross, marker_size);
            ImPlot::PlotScatter("IQ", &buf[0], &buf[1], N, 0, 0, sizeof(std::complex<float>));
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Time Plot")) {
            ImPlot::SetupAxisLinks(ImAxis_Y1, &y_min, &y_max);
            ImPlot::SetupAxis(ImAxis_X1, NULL, ImPlotAxisFlags_AutoFit);
            ImPlot::PlotLine("I", &buf[0], N, 1, 0, 0, 0, sizeof(std::complex<float>));
            ImPlot::PlotLine("Q", &buf[1], N, 1, 0, 0, 0, sizeof(std::complex<float>));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}
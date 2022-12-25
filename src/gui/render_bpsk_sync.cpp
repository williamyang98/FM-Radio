#include "render_bpsk_sync.h"

#include "demod/bpsk_synchroniser.h"

#include <imgui.h>
#include <implot.h>
#include <fmt/core.h>
#include "render_util.h"

void Render_BPSK_Sync(BPSK_Synchroniser& sync) {
    if (ImGui::Begin("PLL Symbols")) {
        PlotConstellationSubplot(sync.GetPLLSymbols(), "###PLL Symbols");
    };
    ImGui::End();

    if (ImGui::Begin("Triggers")) {
        if (ImPlot::BeginPlot("Triggers")) {
            auto range = ImPlotRange(-0.5f, 1.5f);
            auto x0 = sync.GetZeroCrossings();
            auto x1 = sync.GetIntDumpTriggers();
            auto label0 = "Zero Crossing Detector";
            auto label1 = "Integrate & Dump";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotStems(label0, reinterpret_cast<const uint8_t*>(x0.data()), (int)x0.size());
            ImPlot::PlotStems(label1, reinterpret_cast<const uint8_t*>(x1.data()), (int)x1.size());
            ImPlot::EndPlot();
        }
    };
    ImGui::End();

    if (ImGui::Begin("Timing Error Detector")) {
        if (ImPlot::BeginPlot("Phase Error")) {
            auto range = ImPlotRange(-1.5f, 1.5f);
            auto x0 = sync.GetTEDRawPhaseError();
            auto x1 = sync.GetTEDPIPhaseError();
            auto label0 = "Raw";
            auto label1 = "PI controller";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label0, x0.data(), (int)x0.size());
            ImPlot::PlotLine(label1, x1.data(), (int)x1.size());
            ImPlot::EndPlot();
        }
    };
    ImGui::End();

    if (ImGui::Begin("Integrate & Dump Filter")) {
        if (ImPlot::BeginPlot("Integrator & Dump Filter")) {
            auto range = ImPlotRange(-1.5f, 1.5f);
            auto x = sync.GetIntDumpFilter();
            auto label = "Integrate & Dump";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            PlotComplexLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    };
    ImGui::End();

    if (ImGui::Begin("PLL##bpsk_plot")) {
        if (ImPlot::BeginPlot("Phase Error")) {
            auto range = ImPlotRange(-1.5f, 1.5f);
            auto x0 = sync.GetPLLRawPhaseError();
            auto x1 = sync.GetPLLPIPhaseError();
            auto label0 = "Raw";
            auto label1 = "PI controller";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label0, x0.data(), (int)x0.size());
            ImPlot::PlotLine(label1, x1.data(), (int)x1.size());
            ImPlot::EndPlot();
        }
    };
    ImGui::End();
}
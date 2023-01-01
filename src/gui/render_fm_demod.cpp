#include "render_fm_demod.h"

#include "demod/broadcast_fm_demod.h"
#include "dsp/calculate_fft_mag.h"
#include <fmt/core.h>

#include <imgui.h>
#include <implot.h>
#include "render_util.h"
#include "render_bpsk_sync.h"

void Render_FM_Demod_Controls(Broadcast_FM_Demod_Controls& controls);
void Render_FM_Demod_Other_Plots(Broadcast_FM_Demod& demod);
void Render_FM_Demod_Audio_Plots(Broadcast_FM_Demod& demod);
void Render_FM_Demod_RDS_Plots(Broadcast_FM_Demod& demod);
void Render_FM_Demod_Pilot_Plots(Broadcast_FM_Demod& demod);

void Render_FM_Demod_Spectrums(Broadcast_FM_Demod& demod);
void Render_Magnitude_Spectrum_Controls(Calculate_FFT_Mag& calc);
// Custom plot for rendering magnitude spectrum of FM demodulator output signal
void Render_FM_Demod_Magnitude_Spectrum(tcb::span<const float> x, const float Fs, const Broadcast_FM_Demod_Analog_Parameters& config, const ImPlotRange range, const char* label);

// Stereo frame data
static void PlotFrameLine(const char *label, const Frame<float>* x, const int N) {
    auto* x0 = reinterpret_cast<const float*>(x);
    auto L_label = fmt::format("{:s} (L)", label);
    auto R_label = fmt::format("{:s} (R)", label);
    ImPlot::PlotLine(L_label.c_str(), &x0[0], N, 1, 0, 0, 0, sizeof(std::complex<float>));
    ImPlot::PlotLine(R_label.c_str(), &x0[1], N, 1, 0, 0, 0, sizeof(std::complex<float>));
}

void Render_FM_Demod(Broadcast_FM_Demod& demod) {
    if (ImGui::Begin("FM Demodulator Controls")) {
        Render_FM_Demod_Controls(demod.GetControls());
        ImGui::Text("Audio L-R Phase Error: %.2f rads", demod.GetAudioLMRPhaseError());
    }
    ImGui::End();

    if (ImGui::Begin("FM Demodulator Other Signals")) {
        auto id = ImGui::GetID("Dockspace FM Demod Other Signals");
        ImGui::DockSpace(id);
        ImGui::PushID(id);
        Render_FM_Demod_Other_Plots(demod);
        ImGui::PopID();
    }
    ImGui::End();

    if (ImGui::Begin("FM Demodulator Audio Signals")) {
        auto id = ImGui::GetID("Dockspace FM Demod Audio Signals");
        ImGui::DockSpace(id);
        ImGui::PushID(id);
        Render_FM_Demod_Audio_Plots(demod);
        ImGui::PopID();
    }
    ImGui::End();

    if (ImGui::Begin("FM Demodulator RDS Signals")) {
        auto id = ImGui::GetID("Dockspace FM Demod RDS Signals");
        ImGui::DockSpace(id);
        ImGui::PushID(id);
        Render_FM_Demod_RDS_Plots(demod);
        ImGui::PopID();
    }
    ImGui::End();
    
    if (ImGui::Begin("FM Demodulator Pilot Signals")) {
        auto id = ImGui::GetID("Dockspace FM Demod Pilot Signals");
        ImGui::DockSpace(id);
        ImGui::PushID(id);
        Render_FM_Demod_Pilot_Plots(demod);
        ImGui::PopID();
    }
    ImGui::End();

    if (ImGui::Begin("FM Demodulator Spectrums")) {
        auto id = ImGui::GetID("Dockspace FM Demod Spectrums");
        ImGui::DockSpace(id);
        ImGui::PushID(id);
        Render_FM_Demod_Spectrums(demod);
        ImGui::PopID();
    }
    ImGui::End();

    if (ImGui::Begin("RDS BPSK Demodulator")) {
        auto id = ImGui::GetID("Dockspace RDS BPSK");
        ImGui::DockSpace(id);
        ImGui::PushID(id);
        Render_BPSK_Sync(demod.GetBPSKSync());
        ImGui::PopID();
    }
    ImGui::End();
}

void Render_FM_Demod_Spectrums(Broadcast_FM_Demod& demod) {
    if (ImGui::Begin("Baseband Spectrum")) {
        if (ImPlot::BeginPlot("Baseband Spectrum")) {
            auto x = demod.GetBasebandMagnitudeSpectrum();
            const float Fs = (float)demod.GetBasebandSampleRate();
            const int N = (int)x.size();
            auto range = ImPlotRange(-80, 20);
            auto label = "Baseband Spectrum";
            const float xscale = Fs/(float)N;
            const float xstart = -Fs/2.0f;

            const int alpha = 125;
            const auto COL_RED   = ImColor(255,0,0,alpha);
            const auto COL_GREEN = ImColor(0,125,0,alpha);
            //const auto COL_BLUE  = ImColor(0,0,255,alpha);

            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size(), xscale, xstart);
            double xline_0 = 0;
            double xline_1 = +(double)demod.GetFMInSampleRate()/2.0f;
            double xline_2 = -(double)demod.GetFMInSampleRate()/2.0f;
            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_1, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_2, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
        Render_Magnitude_Spectrum_Controls(demod.GetBasebandMagnitudeSpectrumControls());
    }
    ImGui::End();

    if (ImGui::Begin("FM Input Spectrum")) {
        if (ImPlot::BeginPlot("FM Input Spectrum")) {
            auto x = demod.GetFMInMagnitudeSpectrum();
            const float Fs = (float)demod.GetFMInSampleRate();
            const int N = (int)x.size();
            auto range = ImPlotRange(-80, 20);
            auto label = "FM Input";
            const float xscale = Fs/(float)N;
            const float xstart = -Fs/2.0f;

            const int alpha = 125;
            const auto COL_RED   = ImColor(255,0,0,alpha);
            // const auto COL_GREEN = ImColor(0,125,0,alpha);
            // const auto COL_BLUE  = ImColor(0,0,255,alpha);

            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size(), xscale, xstart);
            double xline_0 = 0;

            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
        Render_Magnitude_Spectrum_Controls(demod.GetFMInputMagnitudeSpectrumControls());
    }
    ImGui::End();

    if (ImGui::Begin("FM Output Spectrum")) {
        Render_FM_Demod_Magnitude_Spectrum(
            demod.GetFMOutMagnitudeSpectrum(), 
            (float)demod.GetFMOutSampleRate(), 
            demod.GetAnalogParams(),
            ImPlotRange(-100, -20), 
            "FM Out Spectrum"
        );
        Render_Magnitude_Spectrum_Controls(demod.GetSignalMagnitudeSpectrumControls());
    }
    ImGui::End();

    if (ImGui::Begin("Pilot Tone Spectrum")) {
        if (ImPlot::BeginPlot("Pilot Tone Spectrum")) {
            auto x0 = demod.GetPilotMagnitudeSpectrum();
            auto x1 = demod.GetPLLPilotMagnitudeSpectrum();

            const float Fs = (float)demod.GetFMOutSampleRate();
            const float Fc = (float)demod.GetAnalogParams().F_pilot;

            const int N = (int)x0.size();
            auto range = ImPlotRange(-140, 0);
            auto label0 = "Pilot Tone";
            auto label1 = "PLL Pilot Tone";
            const float xscale = Fs/(float)N;
            const float xstart = -Fs/2.0f;

            const int alpha = 125;
            const auto COL_RED   = ImColor(255,0,0,alpha);
            const auto COL_GREEN = ImColor(0,125,0,alpha);
            // const auto COL_BLUE  = ImColor(0,0,255,alpha);

            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label0, x0.data(), (int)x0.size(), xscale, xstart);
            ImPlot::PlotLine(label1, x1.data(), (int)x1.size(), xscale, xstart);
            double xline_0 = 0;
            double xline_1 = -Fc;
            double xline_2 = +Fc;

            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_1, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_2, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
        Render_Magnitude_Spectrum_Controls(demod.GetPilotMagnitudeSpectrumControls());
        Render_Magnitude_Spectrum_Controls(demod.GetPLLPilotMagnitudeSpectrumControls());
    }
    ImGui::End();

    if (ImGui::Begin("Audio L+R Spectrum")) {
        if (ImPlot::BeginPlot("Audio L+R Spectrum")) {
            auto x = demod.GetAudioLPRMagnitudeSpectrum();
            const float Fs = (float)demod.GetAudioSampleRate();
            const float Fc = (float)demod.GetAnalogParams().F_audio_lpr;
            const int N = (int)x.size();
            auto range = ImPlotRange(-140, 0);
            auto label = "Audio L+R Spectrum";
            const float xscale = Fs/(float)N;
            const float xstart = -Fs/2.0f;

            const int alpha = 125;
            const auto COL_RED   = ImColor(255,0,0,alpha);
            const auto COL_GREEN = ImColor(0,125,0,alpha);
            // const auto COL_BLUE  = ImColor(0,0,255,alpha);

            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size(), xscale, xstart);
            double xline_0 = 0;
            double xline_1 = -Fc;
            double xline_2 = +Fc;

            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_1, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_2, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
        Render_Magnitude_Spectrum_Controls(demod.GetAudioLPRMagnitudeSpectrumControls());
    }
    ImGui::End();

    if (ImGui::Begin("Audio L-R Spectrum")) {
        if (ImPlot::BeginPlot("Audio L-R Spectrum")) {
            auto x = demod.GetAudioLMRMagnitudeSpectrum();
            const float Fs = (float)demod.GetAudioSampleRate();
            const float Fc = (float)demod.GetAnalogParams().F_audio_lmr_bandwidth;
            const int N = (int)x.size();
            auto range = ImPlotRange(-140, 0);
            auto label = "Audio L-R Spectrum";
            const float xscale = Fs/(float)N;
            const float xstart = -Fs/2.0f;

            const int alpha = 125;
            const auto COL_RED   = ImColor(255,0,0,alpha);
            const auto COL_GREEN = ImColor(0,125,0,alpha);
            // const auto COL_BLUE  = ImColor(0,0,255,alpha);

            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size(), xscale, xstart);
            double xline_0 = 0;
            double xline_1 = -Fc;
            double xline_2 = +Fc;

            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_1, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_2, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
        Render_Magnitude_Spectrum_Controls(demod.GetAudioLMRMagnitudeSpectrumControls());
    }
    ImGui::End();

    if (ImGui::Begin("RDS Spectrum")) {
        if (ImPlot::BeginPlot("RDS Spectrum")) {
            auto x = demod.GetRDSMagnitudeSpectrum();
            const float Fs = (float)demod.GetRDSSampleRate();
            const float Fc = (float)demod.GetAnalogParams().F_rds_bandwidth;
            const int N = (int)x.size();
            auto range = ImPlotRange(-140, 0);
            auto label = "RDS Spectrum";
            const float xscale = Fs/(float)N;
            const float xstart = -Fs/2.0f;

            const int alpha = 125;
            const auto COL_RED   = ImColor(255,0,0,alpha);
            const auto COL_GREEN = ImColor(0,125,0,alpha);
            // const auto COL_BLUE  = ImColor(0,0,255,alpha);

            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size(), xscale, xstart);
            double xline_0 = 0;
            double xline_1 = -Fc;
            double xline_2 = +Fc;

            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_1, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_2, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
        Render_Magnitude_Spectrum_Controls(demod.GetRDSMagnitudeSpectrumControls());
    }
    ImGui::End();
}

void Render_FM_Demod_Controls(Broadcast_FM_Demod_Controls& controls) {
    using AudioOut = Broadcast_FM_Demod_Controls::AudioOut;
    static auto GetOutputLabel = [](AudioOut mode) {
        switch (mode) {
        case AudioOut::STEREO:  return "Stereo";
        case AudioOut::LMR:     return "L-R";
        case AudioOut::LPR:     return "L+R";
        default:                return "ERROR: UNKNOWN";
        }
    };

    const auto curr_output_mode = controls.audio_out;
    static auto RenderOutputSelectable = [&controls, &curr_output_mode](AudioOut mode) {
        const bool is_selected = curr_output_mode == mode;
        auto label = GetOutputLabel(mode);
        if (ImGui::Selectable(label, is_selected)) {
            controls.audio_out = mode;
        }
    };

    if (ImGui::BeginCombo("Audio Output", GetOutputLabel(curr_output_mode))) {
        RenderOutputSelectable(AudioOut::STEREO);
        RenderOutputSelectable(AudioOut::LPR);
        RenderOutputSelectable(AudioOut::LMR);
        ImGui::EndCombo();
    }

    ImGui::SliderFloat("L-R gain", &controls.audio_stereo_mix_factor, 0.0f, 5.0f);
    ImGui::Checkbox("L-R Denoiser", &controls.is_audio_lmr_aggressive_lpf);
}

void Render_Magnitude_Spectrum_Controls(Calculate_FFT_Mag& calc) {
    using Mode = Calculate_FFT_Mag::Mode;
    using Trigger = Calculate_FFT_Mag::Trigger;

    static 
    auto GetModeLabel = [](Mode mode) {
        switch (mode) {
        case Mode::NORMAL:   return "Normal";
        case Mode::AVERAGE:  return "Average";
        case Mode::MAX_HOLD: return "Max Hold";
        default:                                return "Unknown";
        }
    };

    static
    auto GetTriggerLabel = [](Trigger trigger) {
        switch (trigger) {
        case Trigger::ALWAYS: return "Always";
        case Trigger::SINGLE: return "Single";
        default:                                 return "Unknown";
        }
    };

    // If we are rendering the spectrum raise the single update flag
    calc.RaiseSingleTrigger();
    const auto curr_mode = calc.GetMode();
    const auto curr_trig = calc.GetTrigger();
    auto* mode_label = GetModeLabel(curr_mode);
    auto* trig_label = GetTriggerLabel(curr_trig);

    static 
    auto RenderMode = [](Mode mode, Mode curr_mode, Calculate_FFT_Mag& calc) {
        auto* label = GetModeLabel(mode);
        const bool is_selected = (mode == curr_mode);
        if (ImGui::Selectable(label, is_selected)) {
            calc.SetMode(mode);
        }
    };

    static 
    auto RenderTrigger = [](Trigger trig, Trigger curr_trig, Calculate_FFT_Mag& calc) {
        auto* label = GetTriggerLabel(trig);
        const bool is_selected = (trig == curr_trig);
        if (ImGui::Selectable(label, is_selected)) {
            calc.SetTrigger(trig);
        }
    };

    ImGui::PushID((void*)(&calc));
    if (ImGui::BeginCombo("Mode", mode_label)) {
        RenderMode(Mode::NORMAL, curr_mode, calc);
        RenderMode(Mode::AVERAGE, curr_mode, calc);
        RenderMode(Mode::MAX_HOLD, curr_mode, calc);
        ImGui::EndCombo();
    }

    const bool is_average = (curr_mode == Mode::AVERAGE);
    if (is_average) {
        auto& v = calc.GetAverageBeta();
        ImGui::SliderFloat("Update Beta", &v, 0.0f, 1.0f);
    }

    if (ImGui::BeginCombo("Trigger", trig_label)) {
        RenderTrigger(Trigger::SINGLE, curr_trig, calc);
        RenderTrigger(Trigger::ALWAYS, curr_trig, calc);
        ImGui::EndCombo();
    }
    ImGui::PopID();
}

void Render_FM_Demod_Magnitude_Spectrum(tcb::span<const float> x, const float Fs, const Broadcast_FM_Demod_Analog_Parameters& config, const ImPlotRange range, const char* label) {
    if (ImPlot::BeginPlot(label)) {
        const int N = (int)x.size();

        const int M = N/2;
        auto x0 = x.subspan(M);

        const float xscale = Fs/(float)N;
        const float xstart = 0.0f;

        ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
        ImPlot::PlotLine(label, x0.data(), (int)x0.size(), xscale, xstart);

        // audio signal
        double Faudio_lpf = (double)config.F_audio_lpr;
        // pilot signal
        double Fpilot = (double)config.F_pilot;
        double Fpilot_bw = (double)config.F_pilot_deviation;
        double Fpilot_lower = Fpilot - Fpilot_bw;
        double Fpilot_upper = Fpilot + Fpilot_bw;
        // stereo signal
        double Fstereo = (double)config.F_audio_lmr_center;
        double Fstereo_bw = (double)config.F_audio_lmr_bandwidth;
        double Fstereo_lower = Fstereo - Fstereo_bw;
        double Fstereo_upper = Fstereo + Fstereo_bw;
        // rds digital signal
        double Frds = (double)config.F_rds_center;
        double Frds_bw = (double)config.F_rds_bandwidth;
        double Frds_lower = Frds - Frds_bw;
        double Frds_upper = Frds + Frds_bw;

        const int alpha = 125;
        const auto COL_RED      = ImColor(255,0,0,alpha);
        const auto COL_GREEN    = ImColor(0,200,0,alpha);
        const auto COL_BLUE     = ImColor(0,0,255,alpha);

        int line_id = 0;
        ImPlot::DragLineX(line_id++, &Faudio_lpf, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);

        ImPlot::DragLineX(line_id++, &Fpilot, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
        // ImPlot::DragLineX(line_id++, &Fpilot_lower, COL_BLUE, 1.0f, ImPlotDragToolFlags_NoInputs);
        // ImPlot::DragLineX(line_id++, &Fpilot_upper, COL_BLUE, 1.0f, ImPlotDragToolFlags_NoInputs);

        ImPlot::DragLineX(line_id++, &Fstereo, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
        ImPlot::DragLineX(line_id++, &Fstereo_lower, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);
        ImPlot::DragLineX(line_id++, &Fstereo_upper, COL_GREEN, 1.0f, ImPlotDragToolFlags_NoInputs);

        ImPlot::DragLineX(line_id++, &Frds, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
        ImPlot::DragLineX(line_id++, &Frds_lower, COL_BLUE, 1.0f, ImPlotDragToolFlags_NoInputs);
        ImPlot::DragLineX(line_id++, &Frds_upper, COL_BLUE, 1.0f, ImPlotDragToolFlags_NoInputs);
        ImPlot::EndPlot();
    }
}

void Render_FM_Demod_Other_Plots(Broadcast_FM_Demod& demod) {
    if (ImGui::Begin("IQ Signal")) {
        if (ImPlot::BeginPlot("IQ Signal")) {
            auto range = ImPlotRange(-1.0f, 1.0f);
            auto x = demod.GetFMOutIQ();
            auto label = "IQ Signal";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            PlotComplexLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();
}

void Render_FM_Demod_Audio_Plots(Broadcast_FM_Demod& demod) {
    if (ImGui::Begin("Audio Output")) {
        if (ImPlot::BeginPlot("Audio Output")) {
            auto range = ImPlotRange(-1.0f, 1.0f);
            auto x = demod.GetAudioOut();
            auto label = "Audio Output";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            PlotFrameLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();

    if (ImGui::Begin("Audio L+R")) {
        if (ImPlot::BeginPlot("Audio L+R")) {
            auto range = ImPlotRange(-1.0f, 1.0f);
            auto x = demod.GetLPRAudioOutput();
            auto label = "Audio L+R";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();

    if (ImGui::Begin("Audio L-R")) {
        if (ImPlot::BeginPlot("Audio L-R")) {
            auto range = ImPlotRange(-1.0f, 1.0f);
            auto x = demod.GetLMRAudioOutput();
            auto label = "Audio L-R";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();
}

void Render_FM_Demod_RDS_Plots(Broadcast_FM_Demod& demod) {
    if (ImGui::Begin("RDS Signal")) {
        PlotConstellationSubplot(demod.GetRDSOutput(), "###RDS Signal");
    }
    ImGui::End();

    if (ImGui::Begin("RDS Raw Output Symbols")) {
        PlotConstellationSubplot(demod.GetRDSRawSymbols(), "###RDS Raw Output Symbols");
    };
    ImGui::End();

    if (ImGui::Begin("RDS Pred Symbols")) {
        if (ImPlot::BeginPlot("RDS Pred Symbols")) {
            auto x = demod.GetRDSPredSymbols();
            // BPSK symbols are normalised to [-1.0f, +1.0f] with noise added
            const float x_width = 2.0f;             
            const float x_range = 2.0f*x_width; 
            const float x_step = 0.1f;
            const int total_bins = (int)(x_range/x_step);

            ImPlot::SetupAxisLimits(ImAxis_X1, -x_width, +x_width, ImPlotCond_Once);
            ImPlot::SetupAxis(ImAxis_Y1, NULL, ImPlotAxisFlags_AutoFit);

            ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.5f);
            ImPlot::PlotHistogram("BPSK counts", x.data(), (int)x.size(), total_bins, 1.0f, ImPlotRange(-x_width, +x_width));
            ImPlot::PopStyleVar();

            const int alpha = 125;
            const auto COL_RED = ImColor(255,0,0,alpha);
            const auto COL_BLUE = ImColor(0,0,255,alpha);
            double xline_0 =  0.0f;
            double xline_1 = -1.0f;
            double xline_2 = +1.0f;
            int line_id = 0;
            ImPlot::DragLineX(line_id++, &xline_0, COL_RED, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_1, COL_BLUE, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::DragLineX(line_id++, &xline_2, COL_BLUE, 1.0f, ImPlotDragToolFlags_NoInputs);
            ImPlot::EndPlot();
        }
    }
    ImGui::End();
}

void Render_FM_Demod_Pilot_Plots(Broadcast_FM_Demod& demod) {
    if (ImGui::Begin("Pilot")) {
        if (ImPlot::BeginPlot("Pilot")) {
            auto range = ImPlotRange(-1.5f, 1.5f);
            auto x = demod.GetPilotOutput();
            auto label = "Pilot";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            PlotComplexLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();

    if (ImGui::Begin("PLL##pilot_plot")) {
        if (ImPlot::BeginPlot("PLL")) {
            auto range = ImPlotRange(-1.5f, 1.5f);
            auto x = demod.GetPLLOutput();
            auto label = "PLL";
            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            PlotComplexLine(label, x.data(), (int)x.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();

    if (ImGui::Begin("PLL Phase Error")) {
        if (ImPlot::BeginPlot("PLL Phase Error")) {
            auto range = ImPlotRange(-1.0f, 1.0f);
            auto x0 = demod.Get_PLL_Raw_Phase_Error_Output();
            auto x1 = demod.Get_PLL_LPF_Phase_Error_Output();
            auto label0 = "Raw";
            auto label1 = "LPF";

            ImPlot::SetupAxes("Sample", "Amplitude");
            ImPlot::SetupAxisLimits(ImAxis_Y1, range.Min, range.Max, ImPlotCond_Once);
            ImPlot::PlotLine(label0, x0.data(), (int)x0.size());
            ImPlot::PlotLine(label1, x1.data(), (int)x1.size());
            ImPlot::EndPlot();
        }
    }
    ImGui::End();
}
// Simple SDR app 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <thread>
#include <memory>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include <GLFW/glfw3.h> 
#include <implot.h>
#include <imgui.h>
#include "gui/imgui_skeleton.h"
#include "gui/font_awesome_definitions.h"
#include "gui/render_app.h"
#include "gui/render_portaudio_controls.h"
#include "gui/render_device_selector.h"

#include "app.h"
#include "audio/portaudio_output.h"
#include "audio/resampled_pcm_player.h"
#include "audio/portaudio_utility.h"
#include "device/device_selector.h"
#include "utility/getopt/getopt.h"

class Renderer: public ImguiSkeleton
{
private:
    App& app;
    PaDeviceList& pa_devices;
    PortAudio_Output& pa_output;
    DeviceSelector& device_selector;
public:
    Renderer(App& _app, PaDeviceList& _pa_devices, PortAudio_Output& _pa_output, DeviceSelector& _device_selector)
    : app(_app), pa_devices(_pa_devices), pa_output(_pa_output), device_selector(_device_selector) {}
    virtual GLFWwindow* Create_GLFW_Window(void) {
        return glfwCreateWindow(
            1280, 720, 
            "Broadcast FM Demodulator", 
            NULL, NULL);
    }
    virtual void AfterImguiContextInit() {
        ImPlot::CreateContext();
        ImguiSkeleton::AfterImguiContextInit();
        auto& io = ImGui::GetIO();
        io.IniFilename =  "imgui_main.ini";
        io.Fonts->AddFontFromFileTTF("res/Roboto-Regular.ttf", 15.0f);
        {
            static const ImWchar icons_ranges[] = { ICON_MIN_FA, ICON_MAX_FA };
            ImFontConfig icons_config;
            icons_config.MergeMode = true;
            icons_config.PixelSnapH = true;
            io.Fonts->AddFontFromFileTTF("res/font_awesome.ttf", 16.0f, &icons_config, icons_ranges);
        }
        ImGuiSetupCustomConfig();
    }

    virtual void Render() {
        if (ImGui::Begin("Our workspace")) {
            ImGui::DockSpace(ImGui::GetID("Our workspace"));

            if (ImGui::Begin("Audio Controls")) {
                RenderPortAudioControls(pa_devices, pa_output);
            }
            ImGui::End();

            if (ImGui::Begin("Device Controls")) {
                RenderDeviceSelector(device_selector);
            } 
            ImGui::End();

            RenderApp(app);
        }
        ImGui::End();
    }
    virtual void AfterShutdown() {
        ImPlot::DestroyContext();
    }
};

void usage() {
    fprintf(stderr, 
        "fm_demod_tuner, Demodulate baseband FM signal from sdr device\n\n"
        "\t[-o output filename (default: None)]\n"
        "\t    If no file is provided then stdout is used\n"
        "\t[-b block size (default: 65536)]\n"
        "\t[-h (show usage)]\n"
    );
}

uint32_t power_ceil(uint32_t x) {
    if (x <= 1) return 1;
    uint32_t power = 2;
    x--;
    while (x >>= 1) power <<= 1;
    return power;
};

int main(int argc, char** argv) {
    uint32_t block_size = 65536;
    const char* wr_filename = NULL;

    int opt; 
    while ((opt = getopt_custom(argc, argv, "o:b:h")) != -1) {
        switch (opt) {
        case 'o':
            wr_filename = optarg;
            break;
        case 'b':
            block_size = (uint32_t)(atof(optarg));
            break;
        case 'h':
        default:
            usage();
            return 0;
        }
    }

    block_size = power_ceil(block_size);
    if (block_size <= 0) {
        fprintf(stderr, "Block size must be positive (%d)\n", block_size); 
        return 1;
    }

    FILE* fp_out = stdout;
    if (wr_filename != NULL) {
        fp_out = fopen(wr_filename, "wb+");
        if (fp_out == NULL) {
            fprintf(stderr, "Failed to open file for writing\n");
            return 1;
        }
    }

    fprintf(stderr, "Using a block size of %u\n", block_size);
#ifdef _WIN32
    _setmode(_fileno(fp_out), _O_BINARY);
#endif

    // Setup audio
    auto pa_handler = ScopedPaHandler();
    PaDeviceList pa_devices;
    PortAudio_Output pa_output;
    std::unique_ptr<Resampled_PCM_Player> pcm_player;
    {
        auto& mixer = pa_output.GetMixer();
        auto buf = mixer.CreateManagedBuffer(4);
        auto Fs = pa_output.GetSampleRate();
        pcm_player = std::make_unique<Resampled_PCM_Player>(buf, Fs);

        #ifdef _WIN32
        const auto target_host_api_index = Pa_HostApiTypeIdToHostApiIndex(PORTAUDIO_TARGET_HOST_API_ID);
        const auto target_device_index = Pa_GetHostApiInfo(target_host_api_index)->defaultOutputDevice;
        pa_output.Open(target_device_index);
        #else
        pa_output.Open(Pa_GetDefaultOutputDevice());
        #endif
    }

    // Setup fm demodulator
    auto app = App(block_size);
    app.OnAudioBlock().Attach([&pcm_player](tcb::span<const Frame<float>> x, int Fs) {
        pcm_player->SetInputSampleRate(Fs);
        pcm_player->ConsumeBuffer(x);
    });

    // Setup input
    auto device_selector = DeviceSelector();
    device_selector.OnDeviceChange().Attach([&](Device* device) {
        if (device == NULL) return;
        device->OnData().Attach([&](tcb::span<const std::complex<uint8_t>> x) {
            app.Process(x);
        });
        // ABC classic is our default frequency
        device->SetCenterFrequency(96'900'000);
        device->SetSamplingFrequency(1'024'000);
    });

    auto init_command_thread = std::thread([&device_selector]() {
		device_selector.SearchDevices();
		if (device_selector.GetDeviceList().size() > 0) {
			device_selector.SelectDevice(0);
		}
	});

    // Setup gui
    auto renderer = Renderer(app, pa_devices, pa_output, device_selector);
    const int rv = RenderImguiSkeleton(&renderer);
    init_command_thread.join();
    return rv;
}
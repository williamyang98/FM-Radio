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

#include "app.h"
#include "audio/portaudio_output.h"
#include "audio/resampled_pcm_player.h"
#include "audio/portaudio_utility.h"
#include "getopt/getopt.h"

class Renderer: public ImguiSkeleton
{
private:
    App& app;
    PaDeviceList& pa_devices;
    PortAudio_Output& pa_output;
public:
    Renderer(App& _app, PaDeviceList& _pa_devices, PortAudio_Output& _pa_output)
    : app(_app), pa_devices(_pa_devices), pa_output(_pa_output) {}
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

            RenderApp(app);
        }
        ImGui::End();
    }
    virtual void AfterShutdown() {
        ImPlot::DestroyContext();
    }
};

struct Arguments {
    const char* input_filename = nullptr;
    uint32_t block_size = 65536;
};

void usage() {
    const auto args = Arguments();

    fprintf(stderr, 
        "fm_demod_no_tuner, Demodulate baseband FM signal from input stream\n\n"
        "\t[-i input filename (default: %s)]\n"
        "\t    If no file is provided then stdin is used\n"
        "\t[-b block size (default: %u)]\n"
        "\t[-h (show usage)]\n",
        args.input_filename ? args.input_filename : "None",
        args.block_size
    );
}

uint32_t power_ceil(uint32_t x) {
    if (x <= 1) return 1;
    uint32_t power = 2;
    x--;
    while (x >>= 1) power <<= 1;
    return power;
};

Arguments parse_args(int argc, char** argv) {
    Arguments args;
    int block_size = int(args.block_size);

    int opt; 
    while ((opt = getopt_custom(argc, argv, "i:b:h")) != -1) {
        switch (opt) {
        case 'i':
            args.input_filename = optarg;
            break;
        case 'b':
            block_size = int(atof(optarg));
            break;
        case 'h':
        default:
            usage();
            exit(0);
        }
    }
    
    if (block_size <= 0) {
        fprintf(stderr, "Block size must be positive (%d)\n", block_size);
        exit(1);
    }

    args.block_size = power_ceil(uint32_t(block_size));
    return args;
}

int main(int argc, char** argv) {
    const auto args = parse_args(argc, argv);

    FILE* fp_in = stdin;
    if (args.input_filename != NULL) {
        fp_in = fopen(args.input_filename, "rb");
        if (fp_in == NULL) {
            fprintf(stderr, "Failed to open file for reading (%s)\n", args.input_filename);
            return 1;
        }
    }

    fprintf(stderr, "Using a block size of %u\n", args.block_size);
#ifdef _WIN32
    _setmode(_fileno(fp_in), _O_BINARY);
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
    auto app = App(args.block_size);
    app.OnAudioBlock().Attach([&pcm_player](tcb::span<const Frame<float>> x, int Fs) {
        pcm_player->SetInputSampleRate(Fs);
        pcm_player->ConsumeBuffer(x);
    });

    // Setup input
    bool is_running = true;
    auto input_buf = std::vector<std::complex<uint8_t>>(args.block_size);
    auto runner_thread = std::make_unique<std::thread>([&]() {
        while (is_running) {
            const size_t nb_read = fread((void*)input_buf.data(), sizeof(std::complex<uint8_t>), input_buf.size(), fp_in);
            if (nb_read != input_buf.size()) {
                fprintf(stderr, "Failed to read %zu/%zu bytes\n", nb_read, input_buf.size());
                break;
            }
            app.Process(input_buf);
        }
        is_running = false;
    });

    // Setup gui
    auto renderer = Renderer(app, pa_devices, pa_output);
    const int rv = RenderImguiSkeleton(&renderer);

    // Close app after gui closes
    is_running = false;
    fclose(fp_in);
    runner_thread->join();

    return rv;
}
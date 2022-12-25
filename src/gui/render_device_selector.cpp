#include "render_device_selector.h"

#include "device/device.h"
#include "device/device_selector.h"

#include <limits.h>
#include <imgui.h>
#include "imgui_extensions.h"
#include <fmt/core.h>

void RenderDeviceControls(Device& device);
bool RenderFrequencyControls(uint32_t& x);

void RenderDeviceSelector(DeviceSelector& app) {
	if (ImGui::Begin("Device Controls")) {
		auto* selected_device = app.GetDevice();
		if (ImGui::BeginTabBar("Device Selector Tabs")) {

			if (ImGui::BeginTabItem("Devices")) {
				if (ImGui::Button("Search")) {
					app.SearchDevices();
				}

				std::string preview_label;
				if (selected_device == NULL) {
					preview_label = "None";
				} else {
					const auto& descriptor = selected_device->GetDescriptor();
					preview_label = fmt::format("[{}] {}",
						descriptor.index, descriptor.product);
				}

				if (ImGui::BeginCombo("Devices", preview_label.c_str())) {
					for (auto& device: app.GetDeviceList()) {
						auto label = fmt::format("[{}] Vendor={} Product={} Serial={}",
							device.index, device.vendor, device.product, device.serial);
						const bool is_selected = (selected_device == NULL) ? 
							false : (selected_device->GetDescriptor().index == device.index);

						if (ImGui::Selectable(label.c_str(), is_selected)) {
							if (is_selected) {
								app.CloseDevice();
							} else {
								app.SelectDevice(device.index);
							}
						}
						if (is_selected) {
							ImGui::SetItemDefaultFocus();
						}
					}
					ImGui::EndCombo();
				}
				ImGui::EndTabItem();
			}

			if (ImGui::BeginTabItem("Controls")) {
				auto lock = std::unique_lock(app.GetDeviceMutex());
				selected_device = app.GetDevice();
				if (selected_device != NULL) {
					RenderDeviceControls(*selected_device);
				}
				ImGui::EndTabItem();
			}

			if (ImGui::BeginTabItem("Errors")) {
				auto lock = std::unique_lock(app.GetDeviceMutex());
				selected_device = app.GetDevice();
				if (selected_device != NULL) {
					auto& errors = selected_device->GetErrorList();
					if (ImGui::BeginListBox("###Errors")) {
						for (auto& error: errors) {
							ImGui::Selectable(error.c_str());
						}
						ImGui::EndListBox();
					}
				}
				ImGui::EndTabItem();
			}

			ImGui::EndTabBar();
		}
	}
	ImGui::End();
}

void RenderDeviceControls(Device& device) {
	std::string preview_label;
	if (!device.GetIsGainManual()) {
		preview_label = "Automatic";
	} else {
		preview_label = fmt::format("{:.1f}dB", device.GetSelectedGain());
	}

	auto& gains = device.GetGainList();
	const auto curr_gain = device.GetSelectedGain();
	int selected_index = -1;

	if (device.GetIsGainManual()) {
		for (int i = 0; i < gains.size(); i++) {
			if (gains[i] == curr_gain) {
				selected_index = i;
				break;
			}
		}
	}

	std::string gain_label = "Automatic";
	if (selected_index >= 0) {
		gain_label = fmt::format("{:.1f}dB", gains[selected_index]);
	}

	if (ImGui::SliderInt("Gain", &selected_index, -1, (int)gains.size()-1, gain_label.c_str())) {
		if (selected_index == -1) {
			device.SetAutoGain();
		} else {
			auto gain = gains[selected_index];
			device.SetGain(gain);
		}
	}

	uint32_t freq = device.GetSelectedFrequency();

	ImGui::BeginGroupPanel("Center Frequency");
	if (RenderFrequencyControls(freq)) {
		device.SetCenterFrequency(freq);
	}
	ImGui::EndGroupPanel();

	ImGui::BeginGroupPanel("Hopping Frequency");
	static uint32_t freq_hop = 800'000;
	RenderFrequencyControls(freq_hop);
	if (ImGui::Button("Hop -")) {
		freq -= freq_hop;
		device.SetCenterFrequency(freq);
	}
	ImGui::SameLine();
	if (ImGui::Button("Hop +")) {
		freq += freq_hop;
		device.SetCenterFrequency(freq);
	}
	ImGui::EndGroupPanel();
}

template <typename T>
T clamp(T x, T min, T max) {
	T y = x;
	y = (y > min) ? y : min;
	y = (y > max) ? max : y;
	return y;
}

bool RenderFrequencyControls(uint32_t& x) {
	// Render up 9GHz 
	constexpr int TOTAL_DIGITS = 10;
	int digits[TOTAL_DIGITS];
	constexpr int limits[TOTAL_DIGITS] = {
		4, 9,9,9, 9,9,9, 9,9,9
	};


	// Generate our digits
	uint32_t y = x;
	uint32_t divisor = 1'000'000'000;
	for (int i = 0; i < 10; i++) {
		const uint32_t multiple = y / divisor;
		y = y - multiple*divisor;
		divisor /= 10;
		digits[i] = (int)multiple;
	}

	const float digit_width = ImGui::CalcTextSize("9").x*1.5f + ImGui::GetStyle().FramePadding.x*2.0f;
	ImGui::PushID((void*)(&x));
	ImGui::PushItemWidth(digit_width);
	ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.0f, ImGui::GetStyle().CellPadding.y));
	bool is_changed = false;
	for (int i = 0; i < TOTAL_DIGITS; i++) {
		const auto flags = ImGuiSliderFlags_AlwaysClamp | ImGuiSliderFlags_ClampOnInput;
		ImGui::PushID(i);
		auto& v = digits[i];

		const bool is_start_group = ((TOTAL_DIGITS-i) % 3) == 0;
		if (is_start_group) {
			ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImGui::GetStyle().CellPadding);
			ImGui::Text(" , ");
			ImGui::PopStyleVar();
			ImGui::SameLine();
		}

		if (ImGui::DragInt("###0", &v, 1.0f, 0, limits[i], "%d", flags)) {
			is_changed = true;
		}

		if (ImGui::IsItemHovered()) {
			const float scroll = ImGui::GetIO().MouseWheel;
			if (scroll > 0.0f) {
				v++;
				is_changed = true;
			}
			if (scroll < 0.0f) {
				v--;
				is_changed = true;
			}
		}

		ImGui::PopID();
		if (i != (TOTAL_DIGITS-1)) {
			ImGui::SameLine();
		}
	}
	ImGui::PopStyleVar();
	ImGui::PopItemWidth();
	ImGui::PopID();

	// Recalculate and update x via reference
	if (is_changed) {
		uint64_t new_x = 0;
		uint64_t multiple = 1'000'000'000;
		for (int i = 0; i < TOTAL_DIGITS; i++) {
			new_x += (uint64_t)digits[i] * multiple;
			multiple /= 10;
		}

		new_x = clamp(new_x, (uint64_t)0, (uint64_t)UINT32_MAX);
		x = (uint32_t)new_x;
	}

	return is_changed;
}
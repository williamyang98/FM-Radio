#include "render_rds_database.h"

#include "rds_decoder/rds_database.h"
#include "rds_decoder/rds_programme_type_names.h"

#include <imgui.h>
#include <fmt/core.h>

void Render_RDS_Database(RDS_Database& db) {
    ImGuiTableFlags flags = ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_Borders;
    if (ImGui::BeginTable("Ensemble description", 2, flags)) {
        ImGui::TableSetupColumn("Field", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        #define FIELD_MACRO(name, fmt, ...) {\
            ImGui::PushID(row_id++);\
            ImGui::TableNextRow();\
            ImGui::TableSetColumnIndex(0);\
            ImGui::TextWrapped(name);\
            ImGui::TableSetColumnIndex(1);\
            ImGui::TextWrapped(fmt, __VA_ARGS__);\
            ImGui::PopID();\
        }\

        int row_id = 0;
        FIELD_MACRO("Programme Identifier", "%04X", db.PI_code);
        FIELD_MACRO("Service Name", "'%.*s'", (int)sizeof(db.service_name), db.service_name);
        FIELD_MACRO("Programme Type", "'%s'", PROGRAMME_TYPES[db.programme_type].display_long);
        FIELD_MACRO("Programme Type Name", "'%.*s'", (int)sizeof(db.programme_type_name), db.programme_type_name);
        FIELD_MACRO("Radio Text", "'%.*s'", (int)sizeof(db.radio_text), db.radio_text);
        FIELD_MACRO("Date", "%02d/%02d/%04d", db.datetime.day, db.datetime.month, db.datetime.year);
        FIELD_MACRO("Time", "%02u:%02u", db.datetime.hour, db.datetime.minute);
        FIELD_MACRO("Local Time Offset", "%d", (int)db.local_time_offset);
        FIELD_MACRO("Audio Type", "%s", db.is_music ? "Music" : "Speech");
        FIELD_MACRO("Audio Channels", "%s", db.is_stereo ? "Stereo" : "Mono");
        FIELD_MACRO("Audio Compressed", "%s", db.is_compressed ? "Yes" : "No");
        FIELD_MACRO("Audio Artificial Head", "%s", db.is_artificial_head ? "Yes" : "No");
        FIELD_MACRO("Dynamic Programme Type", "%s", db.is_dynamic_program_type ? "Yes" : "No");
        #undef FIELD_MACRO

        ImGui::EndTable();
    }

    if (ImGui::Button("Reset Database")) {
        db.Reset();
    }
}
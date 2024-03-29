#include "render_app.h"
#include "app.h"

#include <fmt/core.h>

#include <imgui.h>
#include <implot.h>
#include "render_fm_demod.h"
#include "render_rds_database.h"

void RenderApp(App& app) {
    ImGui::PushID("FM Demod GUI");
    Render_FM_Demod(app.GetFMDemod());
    ImGui::PopID();

    ImGui::PushID("RDS DB");
    if (ImGui::Begin("RDS Database")) {
        Render_RDS_Database(app.GetRDSDatabase());
    }
    ImGui::End();
    ImGui::PopID();
}
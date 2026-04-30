<script setup lang="ts">
import { computed, ref } from "vue";
import { useMutation, useQuery, useQueryClient } from "@tanstack/vue-query";
import { ElMessage, ElMessageBox } from "element-plus";

import { buildApiUrl, cancelRun, deleteRun, getRunLogs, listRuns, packageArtifacts } from "../lib/api";
import type { RunStatus } from "../generated/api";

const queryClient = useQueryClient();
const selectedRun = ref<RunStatus | null>(null);
const logRunId = ref<string | null>(null);
const detailsDialogOpen = ref(false);
const logsDialogOpen = ref(false);

const runsQuery = useQuery({
  queryKey: ["runs"],
  queryFn: listRuns,
  refetchInterval: 5000,
});

const logsQuery = useQuery({
  queryKey: ["run-logs", logRunId],
  queryFn: () => getRunLogs(logRunId.value as string, 300),
  enabled: computed(() => Boolean(logRunId.value)),
  refetchInterval: computed(() => (logRunId.value ? 5000 : false)),
});

const cancelMutation = useMutation({
  mutationFn: (runId: string) => cancelRun(runId),
  onSuccess: async (run) => {
    ElMessage.success(`Cancelled ${run.run_id}`);
    await queryClient.invalidateQueries({ queryKey: ["runs"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to cancel run");
  },
});

const deleteMutation = useMutation({
  mutationFn: (runId: string) => deleteRun(runId),
  onSuccess: async (result) => {
    ElMessage.success(`Deleted ${result.run_id}`);
    if (selectedRun.value?.run_id === result.run_id) {
      selectedRun.value = null;
      detailsDialogOpen.value = false;
    }
    if (logRunId.value === result.run_id) {
      logRunId.value = null;
      logsDialogOpen.value = false;
    }
    await queryClient.invalidateQueries({ queryKey: ["runs"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to delete run");
  },
});

const packageMutation = useMutation({
  mutationFn: (runId: string) => packageArtifacts(runId),
  onSuccess: async (result) => {
    ElMessage.success(`Artifacts packaged for ${result.run_id}`);
    await queryClient.invalidateQueries({ queryKey: ["runs"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to package artifacts");
  },
});

const runs = computed(() => runsQuery.data.value?.runs || []);
const logLines = computed(() => logsQuery.data.value?.tail.join("\n") || "No logs yet.");

function statusTagType(status: string) {
  const normalized = status.toLowerCase();
  if (normalized === "completed") return "success";
  if (normalized === "failed") return "danger";
  if (normalized === "cancelled") return "warning";
  return "info";
}

function formatTimestamp(value?: string | null) {
  if (!value) {
    return "—";
  }
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? value : date.toLocaleString();
}

function openLogs(run: RunStatus) {
  logRunId.value = run.run_id;
  selectedRun.value = run;
  logsDialogOpen.value = true;
}

function openDetails(run: RunStatus) {
  selectedRun.value = run;
  detailsDialogOpen.value = true;
}

async function confirmDelete(run: RunStatus) {
  try {
    await ElMessageBox.confirm(
      `Delete ${run.run_id} and its recorded files?`,
      "Delete run",
      { type: "warning" },
    );
    deleteMutation.mutate(run.run_id);
  } catch {
    // User cancelled the dialog.
  }
}

function downloadArtifact(run: RunStatus) {
  window.open(buildApiUrl(`/runs/${encodeURIComponent(run.run_id)}/artifacts/download`), "_blank");
}

function downloadLog(run: RunStatus) {
  window.open(buildApiUrl(`/runs/${encodeURIComponent(run.run_id)}/logs/download`), "_blank");
}
</script>

<template>
  <section class="page-grid">
    <article class="panel panel--strong">
      <div class="panel__body stack">
        <div>
          <h2 class="panel__title">Tracked runs</h2>
          <p class="panel__copy">
            The UI now uses a proper `GET /api/runs` endpoint, which keeps run discovery and
            polling in the backend instead of the browser.
          </p>
        </div>

        <el-table v-if="runs.length" :data="runs" stripe>
          <el-table-column prop="run_id" label="Run ID" min-width="170" />
          <el-table-column prop="run_name" label="Name" min-width="140" />
          <el-table-column prop="workflow" label="Workflow" width="110" />
          <el-table-column prop="runtime" label="Runtime" width="120" />
          <el-table-column label="Status" width="130">
            <template #default="{ row }">
              <el-tag :type="statusTagType(row.status)">{{ row.status }}</el-tag>
            </template>
          </el-table-column>
          <el-table-column label="Submitted" min-width="180">
            <template #default="{ row }">{{ formatTimestamp(row.submitted_at) }}</template>
          </el-table-column>
          <el-table-column label="Message" min-width="220">
            <template #default="{ row }">{{ row.message || "—" }}</template>
          </el-table-column>
          <el-table-column label="Actions" min-width="340" fixed="right">
            <template #default="{ row }">
              <div class="table-actions">
                <el-button size="small" @click="openDetails(row)">Details</el-button>
                <el-button size="small" @click="openLogs(row)">Logs</el-button>
                <el-button size="small" @click="packageMutation.mutate(row.run_id)">Package</el-button>
                <el-button size="small" @click="downloadArtifact(row)">Download</el-button>
                <el-button size="small" @click="downloadLog(row)">Log file</el-button>
                <el-button size="small" type="warning" @click="cancelMutation.mutate(row.run_id)">
                  Cancel
                </el-button>
                <el-button size="small" type="danger" @click="confirmDelete(row)">Delete</el-button>
              </div>
            </template>
          </el-table-column>
        </el-table>

        <div v-else class="empty-state">
          <p>No runs are currently recorded by the backend.</p>
          <p class="muted">Submit a run from the dashboard to populate this table.</p>
        </div>
      </div>
    </article>

    <el-dialog
      v-model="detailsDialogOpen"
      :title="selectedRun ? `Run details: ${selectedRun.run_id}` : 'Run details'"
      width="760px"
    >
      <template v-if="selectedRun">
        <el-descriptions :column="1" border>
          <el-descriptions-item label="Run name">{{ selectedRun.run_name || "—" }}</el-descriptions-item>
          <el-descriptions-item label="Status">{{ selectedRun.status }}</el-descriptions-item>
          <el-descriptions-item label="Submitted">{{ formatTimestamp(selectedRun.submitted_at) }}</el-descriptions-item>
          <el-descriptions-item label="Input directory">{{ selectedRun.input_dir || "—" }}</el-descriptions-item>
          <el-descriptions-item label="Run directory">{{ selectedRun.run_dir || "—" }}</el-descriptions-item>
          <el-descriptions-item label="Command">
            <pre class="code-panel">{{ selectedRun.command || "No command recorded" }}</pre>
          </el-descriptions-item>
        </el-descriptions>
      </template>
    </el-dialog>

    <el-dialog v-model="logsDialogOpen" :title="logRunId ? `Logs: ${logRunId}` : 'Logs'" width="860px">
      <pre class="code-panel">{{ logLines }}</pre>
    </el-dialog>
  </section>
</template>

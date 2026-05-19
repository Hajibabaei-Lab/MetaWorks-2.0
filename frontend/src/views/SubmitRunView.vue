<script setup lang="ts">
import { computed, reactive, watch } from "vue";
import { useRouter } from "vue-router";
import { useMutation, useQuery, useQueryClient } from "@tanstack/vue-query";
import { ElMessage } from "element-plus";
import { ZodError } from "zod";

import { getHealth, getSettingsPaths, listProfiles, submitRun } from "../lib/api";
import { buildRunSubmission, type SubmitRunFormState } from "../lib/validation";

const router = useRouter();
const queryClient = useQueryClient();

const form = reactive<SubmitRunFormState>({
  profile: "coi",
  workflow: "esv",
  run_name: "",
  runtime: "docker",
  container_image: "docker://metaworks:latest",
  bind_paths_text: "",
  cache_dir: "",
  input_dir: "",
  sample_source: "folder",
  samples_csv: "",
  notes: "",
  cores: "",
  dry_run: false,
  keep_outputs: true,
  config_overrides_text: "{\n  \"pipeline\": {\n    \"name\": \"esv\"\n  }\n}",
});

const healthQuery = useQuery({
  queryKey: ["health"],
  queryFn: getHealth,
  refetchInterval: 15000,
});

const settingsQuery = useQuery({
  queryKey: ["settings"],
  queryFn: getSettingsPaths,
});

const profilesQuery = useQuery({
  queryKey: ["profiles"],
  queryFn: listProfiles,
});

const healthStatus = computed(() => healthQuery.data.value || "");
const settings = computed(() => settingsQuery.data.value || null);
const profiles = computed(() => profilesQuery.data.value?.profiles || []);
const isSubmitting = computed(() => submitMutation.isPending.value);

watch(settingsQuery.data, (settings) => {
  if (!settings) {
    return;
  }

  const allowed = settings.allowed_runtimes as Array<SubmitRunFormState["runtime"]>;
  if (allowed.length > 0 && !allowed.includes(form.runtime)) {
    form.runtime = allowed[0];
  }

  if (!form.container_image || form.container_image === "docker://metaworks:latest") {
    form.container_image = settings.container_image;
  }

  if (!form.input_dir) {
    form.input_dir = `${settings.repo_root}/tests/testing_data`;
  }
}, { immediate: true });

watch(profilesQuery.data, (data) => {
  if (data?.profiles.length && !data.profiles.some((profile) => profile.name === form.profile)) {
    form.profile = data.profiles[0].name;
  }
}, { immediate: true });

const submissionPreview = computed(() => {
  try {
    return JSON.stringify(buildRunSubmission(form), null, 2);
  } catch (error) {
    return String(error instanceof Error ? error.message : error);
  }
});

const allowedRuntimes = computed(() => {
  const values = settings.value?.allowed_runtimes || ["docker", "apptainer"];
  return values.filter(
    (value): value is SubmitRunFormState["runtime"] =>
      value === "conda" || value === "docker" || value === "apptainer",
  );
});

const submitMutation = useMutation({
  mutationFn: () => submitRun(buildRunSubmission(form)),
  onSuccess: async (run) => {
    await queryClient.invalidateQueries({ queryKey: ["runs"] });
    ElMessage.success(`Run ${run.run_id} submitted`);
    router.push({ name: "runs" });
  },
  onError: (error) => {
    if (error instanceof ZodError) {
      ElMessage.error(error.issues[0]?.message || "Please fix the form errors");
      return;
    }

    ElMessage.error(error instanceof Error ? error.message : "Submission failed");
  },
});

const runtimeHelp = computed(() => {
  if (form.runtime === "conda") {
    return "Runs Snakemake directly on the backend host.";
  }

  if (form.runtime === "apptainer") {
    return "Uses an Apptainer or Singularity image and optional cache path.";
  }

  return "Uses Docker with the configured container image and bind paths.";
});
</script>

<template>
  <section class="page-grid">
    <article class="panel panel--strong page-grid__wide">
      <div class="panel__body stack">
        <div>
          <h2 class="panel__title">Submit a pipeline run</h2>
          <p class="panel__copy">
            Launch a run against the backend’s `/api` contract. The form mirrors the API payload so
            it stays maintainable even as the execution engine evolves.
          </p>
          <div class="meta-row">
            <span class="meta-pill">
              API health: {{ healthStatus === "ok" ? "healthy" : "checking" }}
            </span>
            <span class="meta-pill">
              Retention: {{ settings?.retention_policy || "loading" }}
            </span>
            <span class="meta-pill">
              Repo root: {{ settings?.repo_root || "loading" }}
            </span>
          </div>
        </div>

        <el-form label-position="top" class="stack">
          <div class="split-grid">
            <el-form-item label="Profile">
              <el-select v-model="form.profile" filterable>
                <el-option
                  v-for="profile in profiles"
                  :key="profile.name"
                  :label="profile.name"
                  :value="profile.name"
                >
                  <div>
                    <strong>{{ profile.name }}</strong>
                    <div class="muted">{{ profile.description || profile.marker || profile.file }}</div>
                  </div>
                </el-option>
              </el-select>
            </el-form-item>

            <el-form-item label="Workflow">
              <el-segmented
                v-model="form.workflow"
                :options="[
                  { label: 'ESV', value: 'esv' },
                  { label: 'OTU', value: 'otu' },
                ]"
              />
            </el-form-item>
          </div>

          <div class="split-grid">
            <el-form-item label="Run name">
              <el-input v-model="form.run_name" placeholder="COI spring batch" />
            </el-form-item>

            <el-form-item label="Runtime">
              <el-select v-model="form.runtime">
                <el-option
                  v-for="runtime in allowedRuntimes"
                  :key="runtime"
                  :label="runtime"
                  :value="runtime"
                />
              </el-select>
              <div class="muted">{{ runtimeHelp }}</div>
            </el-form-item>
          </div>

          <div class="split-grid">
            <el-form-item label="Input directory">
              <el-input v-model="form.input_dir" placeholder="/data/fastq" />
              <div class="muted">
                Quick start: the bundled sample data lives at
                <code>{{ settings?.repo_root }}/tests/testing_data</code>.
              </div>
            </el-form-item>

            <el-form-item label="Sample source">
              <el-radio-group v-model="form.sample_source">
                <el-radio-button value="folder">Folder</el-radio-button>
                <el-radio-button value="csv">CSV</el-radio-button>
              </el-radio-group>
            </el-form-item>
          </div>

          <el-form-item v-if="form.sample_source === 'csv'" label="Samples CSV">
            <el-input v-model="form.samples_csv" placeholder="/data/samples.csv" />
          </el-form-item>

          <div v-if="form.runtime !== 'conda'" class="split-grid">
            <el-form-item label="Container image">
              <el-input
                v-model="form.container_image"
                placeholder="docker://metaworks:latest"
              />
            </el-form-item>

            <el-form-item label="Cache directory">
              <el-input
                v-model="form.cache_dir"
                placeholder="/workspace/runtime/cache"
              />
            </el-form-item>
          </div>

          <div class="split-grid">
            <el-form-item label="Cores override">
              <el-input v-model="form.cores" placeholder="Leave blank for backend default" />
            </el-form-item>
          </div>

          <el-form-item v-if="form.runtime !== 'conda'" label="Bind paths">
            <el-input
              v-model="form.bind_paths_text"
              type="textarea"
              :rows="4"
              placeholder="/data/inputs:/data/inputs"
            />
          </el-form-item>

          <el-form-item label="Config overrides JSON">
            <el-input
              v-model="form.config_overrides_text"
              type="textarea"
              :rows="10"
              placeholder='{"pipeline":{"name":"esv"}}'
            />
          </el-form-item>

          <el-form-item label="Notes">
            <el-input v-model="form.notes" type="textarea" :rows="3" />
          </el-form-item>

          <div class="split-grid">
            <el-form-item>
              <el-checkbox v-model="form.dry_run">Dry run only</el-checkbox>
            </el-form-item>

            <el-form-item>
              <el-checkbox v-model="form.keep_outputs">Keep outputs until manual cleanup</el-checkbox>
            </el-form-item>
          </div>

          <div class="table-actions">
            <el-button type="primary" :loading="isSubmitting" @click="submitMutation.mutate()">
              Submit run
            </el-button>
            <el-button @click="router.push({ name: 'runs' })">View runs</el-button>
            <el-button @click="form.config_overrides_text = '{}'">Reset overrides</el-button>
          </div>
        </el-form>
      </div>
    </article>

    <aside class="panel page-grid__side">
      <div class="panel__body stack">
        <div>
          <h2 class="panel__title">Payload preview</h2>
          <p class="panel__copy">
            The frontend validates locally with Zod, then sends the same request shape the backend
            accepts over `/api/runs`.
          </p>
        </div>
        <pre class="code-panel">{{ submissionPreview }}</pre>
      </div>
    </aside>
  </section>
</template>

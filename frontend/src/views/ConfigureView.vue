<script setup lang="ts">
import { computed, ref, watch } from "vue";
import { useQuery } from "@tanstack/vue-query";

import { listProfiles } from "../lib/api";
import { useConfigSchema } from "../composables/useConfigSchema";
import ConfigForm from "../components/config/ConfigForm.vue";

const profile = ref("coi");
const workflow = ref<"esv" | "otu">("esv");

const profilesQuery = useQuery({
  queryKey: ["profiles"],
  queryFn: listProfiles,
});

const profiles = computed(() => profilesQuery.data.value?.profiles ?? []);

watch(profilesQuery.data, (data) => {
  if (data?.profiles.length && !data.profiles.some((p) => p.name === profile.value)) {
    profile.value = data.profiles[0].name;
  }
}, { immediate: true });

const { schema, isLoading, error } = useConfigSchema(profile, workflow);
</script>

<template>
  <section class="page-grid">
    <article class="panel panel--strong">
      <div class="panel__body stack">
        <div>
          <h2 class="panel__title">Configure Pipeline</h2>
          <p class="panel__copy">
            Select a marker profile and adjust parameters. Sections collapse by default — expand
            to modify. The live preview renders the merged YAML configuration.
          </p>
        </div>

        <div class="split-grid">
          <el-form-item label="Profile">
            <el-select v-model="profile" filterable>
              <el-option
                v-for="p in profiles"
                :key="p.name"
                :label="p.name"
                :value="p.name"
              >
                <div>
                  <strong>{{ p.name }}</strong>
                  <div class="muted">{{ p.marker || p.description || p.file }}</div>
                </div>
              </el-option>
            </el-select>
          </el-form-item>

          <el-form-item label="Workflow">
            <el-segmented
              v-model="workflow"
              :options="[
                { label: 'ESV', value: 'esv' },
                { label: 'OTU', value: 'otu' },
              ]"
            />
          </el-form-item>
        </div>

        <div v-if="error" class="config-error">
          Failed to load schema: {{ error instanceof Error ? error.message : error }}
        </div>

        <div v-else-if="isLoading" class="empty-state">
          <p>Loading configuration schema...</p>
        </div>

        <ConfigForm v-else-if="schema" :schema="schema" />
      </div>
    </article>
  </section>
</template>

<style scoped>
.config-error {
  padding: 16px 20px;
  border: 1px solid #f56c6c;
  border-radius: 12px;
  background: rgba(245, 108, 108, 0.08);
  color: #b42318;
}
</style>

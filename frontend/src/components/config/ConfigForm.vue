<script setup lang="ts">
import { ref, watch } from "vue";
import { useMutation } from "@tanstack/vue-query";

import type { ConfigSchemaResponse, SectionSchema } from "../../generated/api";
import type { FormValues } from "../../composables/configFormState";
import { buildFormDefaults, buildConfigOverrides, isSectionEnabled } from "../../composables/configFormState";
import { renderConfig } from "../../lib/api";

import ModuleToggles from "./ModuleToggles.vue";
import ConfigSection from "./ConfigSection.vue";
import ConfigPreview from "./ConfigPreview.vue";

const props = defineProps<{
  schema: ConfigSchemaResponse;
}>();

const formValues = ref<FormValues>(buildFormDefaults(props.schema.sections));

watch(
  () => props.schema,
  () => {
    formValues.value = buildFormDefaults(props.schema.sections);
  },
);

function updateSection(sectionKey: string, values: Record<string, unknown>) {
  formValues.value = { ...formValues.value, [sectionKey]: values };
}

const renderMutation = useMutation({
  mutationFn: () =>
    renderConfig({
      profile: props.schema.profile,
      workflow: props.schema.workflow as "esv" | "otu",
      overrides: buildOverrides(),
    }),
});

function buildOverrides(): Record<string, unknown> {
  return buildConfigOverrides(formValues.value, props.schema.sections);
}

function getOrderedSections(): [string, SectionSchema][] {
  const sectionOrder = [
    "pipeline",
    "input",
    "modules",
    "trimming",
    "denoising",
    "classification",
    "pseudogene_filtering",
    "stats",
    "output",
  ];

  return sectionOrder
    .filter((key) => key in props.schema.sections)
    .map((key) => [key, props.schema.sections[key]]);
}

function sectionEnabled(sectionKey: string, section: SectionSchema): boolean {
  if (!section.enabled_by) return true;
  return isSectionEnabled(section.enabled_by, formValues.value);
}
</script>

<template>
  <div class="config-form">
    <div class="config-form__sections">
      <template v-for="[sectionKey, section] in getOrderedSections()" :key="sectionKey">
        <ModuleToggles
          v-if="sectionKey === 'modules'"
          :fields="section.fields"
          :form-values="(formValues.modules as Record<string, unknown>)"
          @update:form-values="(v) => updateSection('modules', v)"
        />

        <ConfigSection
          v-else-if="sectionEnabled(sectionKey, section)"
          :label="section.label"
          :fields="section.fields"
          :form-values="(formValues[sectionKey] as Record<string, unknown>)"
          :default-collapsed="section.collapsed ?? false"
          @update:form-values="(v) => updateSection(sectionKey, v)"
        />
      </template>
    </div>

    <ConfigPreview
      :rendered-yaml="renderMutation.data.value?.merged ?? null"
      :is-loading="renderMutation.isPending.value"
      @render="renderMutation.mutate()"
    />
  </div>
</template>

<style scoped>
.config-form {
  display: grid;
  grid-template-columns: 1fr 380px;
  gap: 24px;
  align-items: start;
}

.config-form__sections {
  display: grid;
  gap: 16px;
}

@media (max-width: 1100px) {
  .config-form {
    grid-template-columns: 1fr;
  }
}
</style>

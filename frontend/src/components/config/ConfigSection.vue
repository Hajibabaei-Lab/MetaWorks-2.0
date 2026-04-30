<script setup lang="ts">
import { ref } from "vue";

import type { FieldSchema } from "../../generated/api";
import { evaluateVisibility, type FormValues } from "../../composables/configFormState";

import FieldRenderer from "./FieldRenderer.vue";

const props = defineProps<{
  label: string;
  fields: Record<string, FieldSchema>;
  formValues: Record<string, unknown>;
  visibleWhen?: string;
  defaultCollapsed?: boolean;
}>();

const emit = defineEmits<{
  "update:formValues": [value: Record<string, unknown>];
}>();

const collapsed = ref(props.defaultCollapsed ?? false);

function onFieldInput(fieldKey: string, value: unknown) {
  emit("update:formValues", { ...props.formValues, [fieldKey]: value });
}

function isVisible(field: FieldSchema): boolean {
  if (!field.visible_when) return true;
  const parentValues: FormValues = {};
  return evaluateVisibility(field.visible_when, { modules: props.formValues, ...parentValues });
}

const visibleFields = () =>
  Object.entries(props.fields).filter(([, f]) => isVisible(f));

function sectionVisible(): boolean {
  if (!props.visibleWhen) return true;
  return evaluateVisibility(props.visibleWhen, { modules: props.formValues });
}
</script>

<template>
  <div v-if="sectionVisible()" class="config-section">
    <button class="config-section__header" @click="collapsed = !collapsed">
      <span class="config-section__label">{{ label }}</span>
      <span class="config-section__toggle">{{ collapsed ? "+" : "−" }}</span>
    </button>

    <div v-show="!collapsed" class="config-section__body">
      <template v-for="[key, field] in visibleFields()" :key="key">
        <FieldRenderer
          v-if="field.type !== 'group'"
          :field-key="key"
          :schema="field"
          :model-value="formValues[key]"
          @update:model-value="(v: unknown) => onFieldInput(key, v)"
        />

        <div v-else class="config-group">
          <h4 class="config-group__title">{{ field.label || key }}</h4>
          <template v-if="field.fields">
            <FieldRenderer
              v-for="[subKey, subField] in Object.entries(field.fields)"
              :key="subKey"
              :field-key="subKey"
              :schema="subField"
              :model-value="(formValues[key] as Record<string, unknown>)?.[subKey]"
              @update:model-value="(v: unknown) => {
                const group = { ...((formValues[key] as Record<string, unknown>) ?? {}) };
                group[subKey] = v;
                onFieldInput(key, group);
              }"
            />
          </template>
        </div>
      </template>
    </div>
  </div>
</template>

<style scoped>
.config-section {
  border: 1px solid var(--border);
  border-radius: 16px;
  overflow: hidden;
  background: var(--panel);
}

.config-section__header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  width: 100%;
  padding: 14px 20px;
  border: none;
  background: rgba(21, 70, 160, 0.04);
  cursor: pointer;
  font-size: 1rem;
  font-weight: 600;
  color: var(--text);
  text-align: left;
}

.config-section__header:hover {
  background: rgba(21, 70, 160, 0.08);
}

.config-section__toggle {
  font-size: 1.2rem;
  color: var(--muted);
  width: 24px;
  text-align: center;
}

.config-section__body {
  padding: 16px 20px;
  display: grid;
  gap: 12px;
}

.config-group {
  padding: 12px 16px;
  border: 1px dashed var(--border);
  border-radius: 12px;
}

.config-group__title {
  margin: 0 0 10px;
  font-size: 0.95rem;
  color: var(--accent);
}
</style>

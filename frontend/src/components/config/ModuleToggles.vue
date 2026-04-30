<script setup lang="ts">
import type { FieldSchema } from "../../generated/api";

const props = defineProps<{
  fields: Record<string, FieldSchema>;
  formValues: Record<string, unknown>;
}>();

const emit = defineEmits<{
  "update:formValues": [value: Record<string, unknown>];
}>();

const toggleFields = () =>
  Object.entries(props.fields).filter(([, f]) => f.type === "boolean");

const engineField = () =>
  Object.entries(props.fields).find(([, f]) => f.type === "select");
</script>

<template>
  <div class="module-toggles">
    <div class="module-toggles__switches">
      <template v-for="[key, field] in toggleFields()" :key="key">
        <label class="module-toggle">
          <el-switch
            :model-value="(formValues[key] as boolean)"
            @update:model-value="(v: boolean) => emit('update:formValues', { ...formValues, [key]: v })"
          />
          <span class="module-toggle__label">{{ field.description || key }}</span>
        </label>
      </template>
    </div>

    <template v-if="engineField()">
      <el-form-item :label="engineField()![1].description || 'Engine'">
        <el-select
          :model-value="formValues[engineField()![0]]"
          @update:model-value="(v: unknown) => emit('update:formValues', { ...formValues, [engineField()![0]]: v })"
        >
          <el-option
            v-for="opt in engineField()![1].options ?? []"
            :key="String(opt)"
            :label="String(opt)"
            :value="opt"
          />
        </el-select>
      </el-form-item>
    </template>
  </div>
</template>

<style scoped>
.module-toggles {
  display: grid;
  gap: 12px;
}

.module-toggles__switches {
  display: flex;
  flex-wrap: wrap;
  gap: 20px;
}

.module-toggle {
  display: flex;
  align-items: center;
  gap: 8px;
  cursor: pointer;
}

.module-toggle__label {
  font-weight: 500;
  text-transform: capitalize;
}
</style>

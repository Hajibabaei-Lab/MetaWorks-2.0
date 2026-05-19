<script setup lang="ts">
import type { FieldSchema } from "../../generated/api";

const props = defineProps<{
  fieldKey: string;
  schema: FieldSchema;
  modelValue: unknown;
}>();

const emit = defineEmits<{
  "update:modelValue": [value: unknown];
}>();

function onInput(value: unknown) {
  emit("update:modelValue", value);
}

function getMin(): number | undefined {
  const c = props.schema.constraints;
  if (!c) return undefined;
  return (c.ge as number) ?? (c.gt as number) ?? undefined;
}

function getMax(): number | undefined {
  const c = props.schema.constraints;
  if (!c) return undefined;
  return (c.le as number) ?? (c.lt as number) ?? undefined;
}

function getStep(): number {
  return props.schema.type === "float" ? 0.01 : 1;
}
</script>

<template>
  <el-form-item :label="schema.description || fieldKey">
    <el-input-number
      v-if="schema.type === 'integer'"
      :model-value="(modelValue as number)"
      :min="getMin()"
      :max="getMax()"
      :step="1"
      controls-position="right"
      @update:model-value="onInput"
    />

    <el-input-number
      v-else-if="schema.type === 'float'"
      :model-value="(modelValue as number)"
      :min="getMin()"
      :max="getMax()"
      :step="getStep()"
      :precision="2"
      controls-position="right"
      @update:model-value="onInput"
    />

    <el-switch
      v-else-if="schema.type === 'boolean'"
      :model-value="(modelValue as boolean)"
      @update:model-value="onInput"
    />

    <el-select
      v-else-if="schema.type === 'select'"
      :model-value="modelValue"
      @update:model-value="onInput"
    >
      <el-option
        v-for="opt in schema.options ?? []"
        :key="String(opt)"
        :label="String(opt)"
        :value="opt"
      />
    </el-select>

    <el-input
      v-else-if="schema.type === 'file_ref'"
      :model-value="String(modelValue ?? '')"
      placeholder="Path to file..."
      @update:model-value="onInput"
    >
      <template #prefix>
        <el-icon><span class="file-icon">&#128196;</span></el-icon>
      </template>
    </el-input>

    <el-input
      v-else
      :model-value="String(modelValue ?? '')"
      @update:model-value="onInput"
    />
  </el-form-item>
</template>

<style scoped>
.file-icon {
  font-style: normal;
}
</style>

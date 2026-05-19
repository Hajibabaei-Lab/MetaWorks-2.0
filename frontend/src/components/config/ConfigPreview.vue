<script setup lang="ts">
defineProps<{
  renderedYaml: string | null;
  isLoading: boolean;
}>();

const emit = defineEmits<{
  render: [];
}>();
</script>

<template>
  <aside class="config-preview">
    <div class="config-preview__header">
      <h3 class="config-preview__title">Live Preview</h3>
      <el-button size="small" type="primary" :loading="isLoading" @click="emit('render')">
        Render YAML
      </el-button>
    </div>
    <pre class="code-panel config-preview__code">{{ renderedYaml ?? "Click Render YAML to preview the merged configuration." }}</pre>
  </aside>
</template>

<style scoped>
.config-preview {
  position: sticky;
  top: 24px;
  border: 1px solid var(--border);
  border-radius: 16px;
  background: var(--panel);
  overflow: hidden;
}

.config-preview__header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 14px 20px;
  background: rgba(21, 70, 160, 0.04);
}

.config-preview__title {
  margin: 0;
  font-size: 1rem;
  font-weight: 600;
}

.config-preview__code {
  min-height: 400px;
  max-height: 70vh;
  overflow: auto;
  margin: 0;
  border-radius: 0;
}
</style>

<script setup lang="ts">
import { computed } from "vue";
import { useMutation, useQuery, useQueryClient } from "@tanstack/vue-query";
import { ElMessage, ElMessageBox } from "element-plus";

import {
  deleteAdapter,
  deleteClassifier,
  listAdapters,
  listClassifiers,
  uploadAdapter,
  uploadClassifier,
} from "../lib/api";

const queryClient = useQueryClient();

const adaptersQuery = useQuery({
  queryKey: ["adapters"],
  queryFn: listAdapters,
});

const classifiersQuery = useQuery({
  queryKey: ["classifiers"],
  queryFn: listClassifiers,
});

const uploadAdapterMutation = useMutation({
  mutationFn: uploadAdapter,
  onSuccess: async (result) => {
    ElMessage.success(`Uploaded ${result.name}`);
    await queryClient.invalidateQueries({ queryKey: ["adapters"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to upload adapter");
  },
});

const uploadClassifierMutation = useMutation({
  mutationFn: uploadClassifier,
  onSuccess: async (result) => {
    ElMessage.success(`Uploaded ${result.name}`);
    await queryClient.invalidateQueries({ queryKey: ["classifiers"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to upload classifier");
  },
});

const deleteAdapterMutation = useMutation({
  mutationFn: deleteAdapter,
  onSuccess: async () => {
    ElMessage.success("Adapter deleted");
    await queryClient.invalidateQueries({ queryKey: ["adapters"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to delete adapter");
  },
});

const deleteClassifierMutation = useMutation({
  mutationFn: deleteClassifier,
  onSuccess: async () => {
    ElMessage.success("Classifier deleted");
    await queryClient.invalidateQueries({ queryKey: ["classifiers"] });
  },
  onError: (error) => {
    ElMessage.error(error instanceof Error ? error.message : "Failed to delete classifier");
  },
});

const adapters = computed(() => adaptersQuery.data.value?.items || []);
const classifiers = computed(() => classifiersQuery.data.value?.items || []);

function onFilePicked(
  event: Event,
  kind: "adapter" | "classifier",
) {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (!file) {
    return;
  }

  if (kind === "adapter") {
    uploadAdapterMutation.mutate(file);
  } else {
    uploadClassifierMutation.mutate(file);
  }

  input.value = "";
}

async function confirmDelete(kind: "adapter" | "classifier", name: string) {
  try {
    await ElMessageBox.confirm(`Delete ${name}?`, "Delete asset", { type: "warning" });
    if (kind === "adapter") {
      deleteAdapterMutation.mutate(name);
    } else {
      deleteClassifierMutation.mutate(name);
    }
  } catch {
    // User cancelled the dialog.
  }
}
</script>

<template>
  <section class="page-grid">
    <article class="panel panel--strong">
      <div class="panel__body stack">
        <div>
          <h2 class="panel__title">Managed assets</h2>
          <p class="panel__copy">
            Upload and remove adapters or classifiers through the API. The backend remains the
            source of truth for where these files live on disk.
          </p>
        </div>

        <div class="split-grid">
          <section class="panel">
            <div class="panel__body stack">
              <div>
                <h3 class="panel__title">Adapters</h3>
                <p class="panel__copy">Upload FASTA adapter files for trimming workflows.</p>
              </div>
              <input type="file" @change="(event) => onFilePicked(event, 'adapter')" />
              <el-table :data="adapters" stripe>
                <el-table-column label="Filename" min-width="240">
                  <template #default="{ row }">{{ row }}</template>
                </el-table-column>
                <el-table-column label="Actions" width="140">
                  <template #default="{ row }">
                    <el-button size="small" type="danger" @click="confirmDelete('adapter', row)">
                      Delete
                    </el-button>
                  </template>
                </el-table-column>
              </el-table>
            </div>
          </section>

          <section class="panel">
            <div class="panel__body stack">
              <div>
                <h3 class="panel__title">Classifiers</h3>
                <p class="panel__copy">Manage RDP or custom classifier assets.</p>
              </div>
              <input type="file" @change="(event) => onFilePicked(event, 'classifier')" />
              <el-table :data="classifiers" stripe>
                <el-table-column label="Filename" min-width="240">
                  <template #default="{ row }">{{ row }}</template>
                </el-table-column>
                <el-table-column label="Actions" width="140">
                  <template #default="{ row }">
                    <el-button size="small" type="danger" @click="confirmDelete('classifier', row)">
                      Delete
                    </el-button>
                  </template>
                </el-table-column>
              </el-table>
            </div>
          </section>
        </div>
      </div>
    </article>
  </section>
</template>

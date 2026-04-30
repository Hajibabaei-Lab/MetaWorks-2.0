import { computed, type Ref } from "vue";
import { useQuery } from "@tanstack/vue-query";

import { getConfigSchema } from "../lib/api";

export function useConfigSchema(profile: Ref<string>, workflow: Ref<string>) {
  const query = useQuery({
    queryKey: computed(() => ["config-schema", profile.value, workflow.value]),
    queryFn: () => getConfigSchema(profile.value, workflow.value),
    enabled: computed(() => Boolean(profile.value)),
  });

  return {
    schema: computed(() => query.data.value ?? null),
    isLoading: query.isPending,
    error: query.error,
    refetch: query.refetch,
  };
}

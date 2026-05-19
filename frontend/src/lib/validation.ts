import { z } from "zod";

import type {
  ConfigSchemaResponse,
  ConfigSectionsResponse,
  ListAssetsResponse,
  PathsResponse,
  ProfileListResponse,
  RenderConfigResponse,
  RunListResponse,
  RunStatus,
  RunSubmissionRequest,
  UploadResponse,
} from "../generated/api";

export const workflowSchema = z.enum(["esv", "otu"]);
export const runtimeSchema = z.enum(["conda", "docker", "apptainer"]);

export const runStatusSchema: z.ZodType<RunStatus> = z.object({
  run_id: z.string(),
  workflow: workflowSchema,
  runtime: runtimeSchema,
  status: z.string(),
  run_name: z.string().nullable().optional(),
  sample_source: z.string().nullable().optional(),
  input_dir: z.string().nullable().optional(),
  samples_csv: z.string().nullable().optional(),
  notes: z.string().nullable().optional(),
  config_overrides: z.record(z.string(), z.unknown()).nullable().optional(),
  scheduler_job_id: z.string().nullable().optional(),
  submitted_at: z.string().nullable().optional(),
  started_at: z.string().nullable().optional(),
  ended_at: z.string().nullable().optional(),
  pid: z.number().nullable().optional(),
  return_code: z.number().nullable().optional(),
  message: z.string().nullable().optional(),
  run_dir: z.string().nullable().optional(),
  log_path: z.string().nullable().optional(),
  command: z.string().nullable().optional(),
  artifact_path: z.string().nullable().optional(),
  keep_outputs: z.boolean().nullable().optional(),
});

export const runListSchema: z.ZodType<RunListResponse> = z.object({
  runs: z.array(runStatusSchema),
});

export const profileListSchema: z.ZodType<ProfileListResponse> = z.object({
  profiles: z.array(
    z.object({
      name: z.string(),
      description: z.string().optional().default(""),
      marker: z.string().optional().default(""),
      file: z.string(),
    }),
  ),
});

export const renderConfigResponseSchema: z.ZodType<RenderConfigResponse> = z.object({
  workflow: workflowSchema,
  merged: z.string(),
});

export const configSectionsResponseSchema: z.ZodType<ConfigSectionsResponse> = z.object({
  workflow: workflowSchema,
  sections: z.record(z.string(), z.unknown()),
});

export const assetsListSchema: z.ZodType<ListAssetsResponse> = z.object({
  items: z.array(z.string()),
});

export const uploadResponseSchema: z.ZodType<UploadResponse> = z.object({
  name: z.string(),
  path: z.string(),
});

export const pathsResponseSchema: z.ZodType<PathsResponse> = z.object({
  repo_root: z.string(),
  runtime_cache: z.string(),
  allowed_runtimes: z.array(z.string()),
  retention_policy: z.string(),
  default_runtime: z.string(),
  container_image: z.string(),
});

const fieldSchema: z.ZodType<import("../generated/api").FieldSchema> = z.lazy(() =>
  z.object({
    type: z.string(),
    default: z.unknown().optional(),
    description: z.string().optional(),
    constraints: z.record(z.string(), z.unknown()).optional(),
    options: z.array(z.unknown()).optional(),
    visible_when: z.string().optional(),
    label: z.string().optional(),
    fields: z.record(z.string(), fieldSchema).optional(),
  }),
);

export const configSchemaResponseSchema: z.ZodType<ConfigSchemaResponse> = z.object({
  profile: z.string(),
  marker: z.string(),
  workflow: z.string(),
  sections: z.record(
    z.string(),
    z.object({
      label: z.string(),
      enabled_by: z.string().optional(),
      collapsed: z.boolean().optional(),
      fields: z.record(z.string(), fieldSchema),
    }),
  ),
});

export type SubmitRunFormState = {
  profile: string;
  workflow: "esv" | "otu";
  run_name: string;
  runtime: "conda" | "docker" | "apptainer";
  container_image: string;
  bind_paths_text: string;
  cache_dir: string;
  input_dir: string;
  sample_source: "folder" | "csv";
  samples_csv: string;
  notes: string;
  cores: string;
  dry_run: boolean;
  keep_outputs: boolean;
  config_overrides_text: string;
};

export const submitRunFormSchema = z
  .object({
    profile: z.string().trim().min(1, "Profile is required"),
    workflow: workflowSchema,
    run_name: z.string().trim().min(1, "Run name is required"),
    runtime: runtimeSchema,
    container_image: z.string().trim().optional(),
    bind_paths_text: z.string(),
    cache_dir: z.string(),
    input_dir: z.string().trim().min(1, "Input directory is required"),
    sample_source: z.enum(["folder", "csv"]),
    samples_csv: z.string(),
    notes: z.string(),
    cores: z.string(),
    dry_run: z.boolean(),
    keep_outputs: z.boolean(),
    config_overrides_text: z.string(),
  })
  .superRefine((value, ctx) => {
    if (
      (value.runtime === "docker" || value.runtime === "apptainer") &&
      !value.container_image?.trim()
    ) {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        message: "Container image is required for Docker and Apptainer runs",
        path: ["container_image"],
      });
    }

    if (value.sample_source === "csv" && !value.samples_csv.trim()) {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        message: "Samples CSV is required when sample source is CSV",
        path: ["samples_csv"],
      });
    }
  });

export function parseJsonObject(text: string, label: string): Record<string, unknown> {
  const raw = text.trim();
  if (!raw) {
    return {};
  }

  let value: unknown;
  try {
    value = JSON.parse(raw);
  } catch (error) {
    throw new Error(`${label} must be valid JSON`);
  }

  if (!value || Array.isArray(value) || typeof value !== "object") {
    throw new Error(`${label} must be a JSON object`);
  }

  return value as Record<string, unknown>;
}

function parseOptionalInt(value: string, label: string): number | null {
  const trimmed = value.trim();
  if (!trimmed) {
    return null;
  }

  const parsed = Number(trimmed);
  if (!Number.isInteger(parsed) || parsed < 1) {
    throw new Error(`${label} must be a positive whole number`);
  }

  return parsed;
}

function parseBindPaths(bindPathsText: string): string[] {
  return bindPathsText
    .split("\n")
    .map((item) => item.trim())
    .filter(Boolean);
}

export function buildRunSubmission(form: SubmitRunFormState): RunSubmissionRequest {
  const parsed = submitRunFormSchema.parse(form);
  const overrides = parseJsonObject(parsed.config_overrides_text, "Config overrides");
  const cores = parseOptionalInt(parsed.cores, "Cores");

  return {
    profile: parsed.profile,
    workflow: parsed.workflow,
    run_name: parsed.run_name,
    runtime: parsed.runtime,
    container_image:
      parsed.runtime === "conda" ? null : (parsed.container_image?.trim() || null),
    bind_paths: parseBindPaths(parsed.bind_paths_text),
    cache_dir: parsed.cache_dir.trim() || null,
    config_overrides: overrides,
    input_dir: parsed.input_dir,
    samples_csv: parsed.sample_source === "csv" ? parsed.samples_csv.trim() : null,
    sample_source: parsed.sample_source,
    notes: parsed.notes.trim() || null,
    cores,
    dry_run: parsed.dry_run,
    keep_outputs: parsed.keep_outputs,
  };
}


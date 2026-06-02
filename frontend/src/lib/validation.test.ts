import { buildRunSubmission, configSchemaResponseSchema, parseJsonObject } from "./validation";

describe("buildRunSubmission", () => {
  it("normalizes the form into an API payload", () => {
    const payload = buildRunSubmission({
      profile: "coi",
      workflow: "esv",
      run_name: "demo",
      runtime: "docker",
      container_image: "docker://metaworks:latest",
      bind_paths_text: "/data:/data\n/work:/work",
      cache_dir: "/tmp/cache",
      input_dir: "/data/fastq",
      sample_source: "folder",
      samples_csv: "",
      notes: "hello",
      cores: "8",
      dry_run: false,
      keep_outputs: true,
      config_overrides_text: "{\"pipeline\":{\"name\":\"esv\"}}",
    });

    expect(payload.bind_paths).toEqual(["/data:/data", "/work:/work"]);
    expect(payload.cores).toBe(8);
    expect(payload.config_overrides).toEqual({ pipeline: { name: "esv" } });
  });

  it("rejects invalid JSON objects", () => {
    expect(() => parseJsonObject("[]", "Overrides")).toThrow("Overrides must be a JSON object");
  });
});

describe("configSchemaResponseSchema", () => {
  it("accepts null optional schema metadata returned by the backend", () => {
    const parsed = configSchemaResponseSchema.parse({
      profile: "coi",
      marker: "COI",
      workflow: "esv",
      sections: {
        pipeline: {
          label: "Pipeline Settings",
          enabled_by: null,
          collapsed: null,
          fields: {
            parallel_jobs: {
              type: "integer",
              default: 4,
              description: null,
              constraints: { ge: 1, le: 32 },
              options: null,
              visible_when: null,
              label: null,
              fields: null,
            },
          },
        },
      },
    });

    expect(parsed.sections.pipeline.enabled_by).toBeUndefined();
    expect(parsed.sections.pipeline.collapsed).toBeUndefined();
    expect(parsed.sections.pipeline.fields.parallel_jobs.description).toBeUndefined();
    expect(parsed.sections.pipeline.fields.parallel_jobs.fields).toBeUndefined();
  });
});

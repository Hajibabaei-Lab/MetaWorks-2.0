import { buildRunSubmission, parseJsonObject } from "./validation";

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
      mem_gb: "32",
      dry_run: false,
      keep_outputs: true,
      config_overrides_text: "{\"pipeline\":{\"name\":\"esv\"}}",
    });

    expect(payload.bind_paths).toEqual(["/data:/data", "/work:/work"]);
    expect(payload.cores).toBe(8);
    expect(payload.mem_gb).toBe(32);
    expect(payload.config_overrides).toEqual({ pipeline: { name: "esv" } });
  });

  it("rejects invalid JSON objects", () => {
    expect(() => parseJsonObject("[]", "Overrides")).toThrow("Overrides must be a JSON object");
  });
});


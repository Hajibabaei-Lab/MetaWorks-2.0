import { buildApiUrl } from "./api";

describe("buildApiUrl", () => {
  it("uses the configured API base", () => {
    expect(buildApiUrl("/runs")).toBe("/api/runs");
    expect(buildApiUrl("health")).toBe("/api/health");
  });
});


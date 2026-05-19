#!/usr/bin/env node

import { execFileSync } from "node:child_process";

const schemaUrl = process.env.SCHEMA_URL || "http://127.0.0.1:8000/api/openapi.json";

execFileSync(
  "npx",
  ["openapi-typescript", schemaUrl, "-o", "src/generated/api.ts"],
  { stdio: "inherit" },
);


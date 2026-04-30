import type { FieldSchema, SectionSchema } from "../generated/api";

export type FormValues = Record<string, unknown>;

export function buildFormDefaults(
  sections: Record<string, SectionSchema>,
): FormValues {
  const result: FormValues = {};

  for (const [sectionKey, section] of Object.entries(sections)) {
    result[sectionKey] = extractSectionDefaults(section.fields);
  }

  return result;
}

function extractSectionDefaults(
  fields: Record<string, FieldSchema>,
): Record<string, unknown> {
  const result: Record<string, unknown> = {};

  for (const [key, field] of Object.entries(fields)) {
    if (field.type === "group" && field.fields) {
      result[key] = extractSectionDefaults(field.fields);
    } else if (field.default !== undefined) {
      result[key] = field.default;
    } else {
      result[key] = null;
    }
  }

  return result;
}

export function buildConfigOverrides(
  formValues: FormValues,
  sections: Record<string, SectionSchema>,
): Record<string, unknown> {
  const result: Record<string, unknown> = {};

  for (const [sectionKey, section] of Object.entries(sections)) {
    const sectionValues = formValues[sectionKey];
    if (!sectionValues || typeof sectionValues !== "object") continue;

    const overrides = diffAgainstDefaults(
      sectionValues as Record<string, unknown>,
      section.fields,
    );

    if (Object.keys(overrides).length > 0) {
      result[sectionKey] = overrides;
    }
  }

  return result;
}

function diffAgainstDefaults(
  values: Record<string, unknown>,
  fields: Record<string, FieldSchema>,
): Record<string, unknown> {
  const result: Record<string, unknown> = {};

  for (const [key, field] of Object.entries(fields)) {
    const value = values[key];

    if (field.type === "group" && field.fields && value && typeof value === "object") {
      const nested = diffAgainstDefaults(
        value as Record<string, unknown>,
        field.fields,
      );
      if (Object.keys(nested).length > 0) {
        result[key] = nested;
      }
      continue;
    }

    const schemaDefault = field.default;
    if (value !== schemaDefault) {
      result[key] = value;
    }
  }

  return result;
}

export function evaluateVisibility(
  visibleWhen: string,
  formValues: FormValues,
): boolean {
  const match = visibleWhen.match(
    /^([a-zA-Z_]+)\.([a-zA-Z_]+)\s*==\s*'(.+)'$/,
  );
  if (!match) return true;

  const [, sectionKey, fieldKey, expectedValue] = match;
  const section = formValues[sectionKey];
  if (!section || typeof section !== "object") return false;

  const actualValue = (section as Record<string, unknown>)[fieldKey];
  return String(actualValue) === expectedValue;
}

export function isSectionEnabled(
  enabledBy: string | undefined,
  formValues: FormValues,
): boolean {
  if (!enabledBy) return true;

  const [sectionKey, fieldKey] = enabledBy.split(".");
  const section = formValues[sectionKey];
  if (!section || typeof section !== "object") return false;

  return Boolean((section as Record<string, unknown>)[fieldKey]);
}

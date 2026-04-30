import { expect, test } from "@playwright/test";

test("renders the shell and primary navigation", async ({ page }) => {
  await page.goto("/");

  await expect(page.getByText("Independent web UI for the MetaWorks runner")).toBeVisible();
  await expect(page.getByRole("link", { name: "Submit Run" })).toBeVisible();
  await expect(page.getByRole("link", { name: "Runs" })).toBeVisible();
  await expect(page.getByRole("link", { name: "Configure" })).toBeVisible();
  await expect(page.getByRole("link", { name: "Assets" })).toBeVisible();
});


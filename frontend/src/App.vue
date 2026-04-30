<script setup lang="ts">
import { RouterLink, RouterView, useRoute } from "vue-router";

const route = useRoute();

const links = [
  { to: "/", label: "Submit Run" },
  { to: "/runs", label: "Runs" },
  { to: "/configure", label: "Configure" },
  { to: "/assets", label: "Assets" },
];
</script>

<template>
  <div class="app-shell">
    <header class="app-header">
      <div class="app-header__content">
        <div>
          <div class="app-kicker">MetaWorks Control Surface</div>
          <h1>Independent web UI for the MetaWorks runner</h1>
          <p>
            The frontend is now a standalone app. The backend stays independently usable through
            its `/api` contract and existing headless workflows.
          </p>
        </div>
        <div class="app-orb" aria-hidden="true"></div>
      </div>
      <nav class="app-nav" aria-label="Primary">
        <RouterLink
          v-for="link in links"
          :key="link.to"
          :to="link.to"
          class="app-nav__link"
          :class="{ 'is-active': route.path === link.to }"
        >
          {{ link.label }}
        </RouterLink>
      </nav>
    </header>

    <main class="app-main">
      <RouterView />
    </main>
  </div>
</template>


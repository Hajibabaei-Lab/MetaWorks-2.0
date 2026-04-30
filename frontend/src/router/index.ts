import { createRouter, createWebHistory } from "vue-router";

import SubmitRunView from "../views/SubmitRunView.vue";
import RunsView from "../views/RunsView.vue";
import ConfigureView from "../views/ConfigureView.vue";
import AssetsView from "../views/AssetsView.vue";

const router = createRouter({
  history: createWebHistory(),
  routes: [
    { path: "/", name: "submit", component: SubmitRunView },
    { path: "/runs", name: "runs", component: RunsView },
    { path: "/configure", name: "configure", component: ConfigureView },
    { path: "/assets", name: "assets", component: AssetsView },
  ],
});

export default router;


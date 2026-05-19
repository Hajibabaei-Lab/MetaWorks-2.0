import { createApp } from "vue";
import { VueQueryPlugin, QueryClient } from "@tanstack/vue-query";
import ElementPlus from "element-plus";
import "element-plus/dist/index.css";

import App from "./App.vue";
import router from "./router";
import "./styles.css";

const app = createApp(App);
const queryClient = new QueryClient();

app.use(router);
app.use(ElementPlus);
app.use(VueQueryPlugin, { queryClient });
app.mount("#app");


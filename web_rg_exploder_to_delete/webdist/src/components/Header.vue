<template>
  <header class="header">
    <h1 class="header__title">
      {{ title }}
    </h1>
    <img
      class="header__logo"
      :src="repliconLogo"
      alt="replicon logo"
    >
  </header>
</template>

<script>
export default {
  props: {
    title: {
      type: String,
      required: true
    }
  },
  data() {
    return {
      repliconLogo: require('../assets/replicon_logo.svg')
    }
  }
}
</script>

<style scoped>
.header {
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
  height: 5em;
  width: 100%;
  margin-bottom: 2rem;
  box-shadow: var(--primary-shadow);
}

.header__title {
  margin-right: 2em;
  color: var(--primary-color);
}

.header__logo {
  height: 3.5em;
}
</style>

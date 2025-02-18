<template>
  <div class="copyright">
    <h3 class="copyright__title">
      {{ copytitle }}
    </h3>
  </div>
</template>

<script>
export default {
  props: {
    journal: {
      type: String,
      default: ''
    },
    readme: {
      type: String,
      default: ''
    },
    copytitle: {
      type: String,
      default: ''
    }
  }
}
</script>

<style scoped>
.copyright {
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
  height: 5em;
  width: 80%;
  margin: 2rem auto;
  border-radius: 8px;
  box-shadow: var(--primary-shadow);
}

.copyright__title {
  margin-right: 2em;
  color: var(--primary-color);
}

</style>

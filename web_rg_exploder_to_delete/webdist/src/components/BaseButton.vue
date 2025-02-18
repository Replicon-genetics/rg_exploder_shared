<template>
  <button
    class="button"
    :disabled="running"
    :label="label"
  >
    {{ running ? '...' : label }}
  </button>
</template>

<script>
export default {
  props: {
    running: {
      type: Boolean,
      default: false
    },
    label: {
      type: String,
      default: ''
    },
  }
}
</script>

<style scoped>
.button {
  background-color: var(--primary-color);
  color: black;
  margin: 0.5rem;
  height: 1.5em;
  width: 10em;
  font-size: 1.0rem;
  font-weight: bold;
  border: none;
  cursor: pointer;
  border-radius: 10%;
  box-shadow: var(--primary-shadow);
  transition: .2s;
  right: 1em;
  bottom: 1em;
}

.button:hover {
  filter: brightness(1.1);
}
</style>

<template>
  <button
    class="button"
    :disabled="running"
  >
    {{ running ? '...' : 'Run' }}
  </button>
</template>

<script>
export default {
  props: {
    running: {
      type: Boolean,
      default: false
    }
  }
}
</script>

<style scoped>
.button {
  background-color: var(--primary-color);
  color: white;
  margin: 0.5rem;
  height: 3em;
  width: 3em;
  font-size: 1.5rem;
  font-weight: bold;
  border: none;
  cursor: pointer;
  border-radius: 50%;
  box-shadow: var(--primary-shadow);
  transition: .2s;
  position: fixed;
  right: 1em;
  bottom: 1em;
}

.button:hover {
  filter: brightness(1.1);
}
</style>

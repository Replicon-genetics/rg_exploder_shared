<template>
  <div class="message-panel">
    <h2>Message panel</h2>

    <h3 v-if="journal">
      Journal
    </h3>
    <pre v-if="journal"> {{ journal }} </pre>

    <h3 v-if="readme">
      Readme
    </h3>
    <pre v-if="readme"> {{ readme }} </pre>
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
    }
  }
}
</script>

<style scoped>
.message-panel {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  width: 80%;
  margin: 2rem auto;
  border-radius: 8px;
  box-shadow: var(--primary-shadow);
}

pre {
  width: 90%;
  white-space: pre-line;
  overflow: scroll;
  overflow-y: scroll;
  border: 1px solid black;
  border-radius: 8px;
  padding: 1em;
  color: black;
  background: white;
}
</style>

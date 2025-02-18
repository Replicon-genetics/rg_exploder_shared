<template>
  <div
    v-if="shouldShow"
    class="checkbox"
  >
    <label v-if="label"> {{ label }} </label>
    <input
      value="modelValue"
      type="checkbox"
      :checked="modelValue"
      :disabled="disabled"
      @change="$emit('update:modelValue', $event.target.checked)"
    >
  </div>
</template>

<script>
export default {
  props: {
    modelValue: {
      type: [Boolean, String],
      default: ''
    },
    label: {
      type: String,
      default: ''
    },
    disabled: {
      type: Boolean,
      default: false
    }
  },
  emits: ['update:modelValue'],
  computed: {
    shouldShow() {
      return (typeof this.modelValue) === 'boolean'
    }
  },
  watch: {
    modelValue: function() {
      console.log(`${this.label} is now ${this.modelValue}`)
    }
  }
}
</script>

<style scoped>
.checkbox {
  max-width: 90%;
  font-size: var(--fnt-normal);
  display: flex;
  justify-content: space-between;
  margin-bottom: 0.5em;
}

input {
  width: 1em;
}
</style>

<template>
  <div class="baseinput">
    <label> {{ label }} </label>
    <div class="input">
      <span> {{ modelValue }} </span>
      <input
        type="range"
        min="0"
        max="100"
        :value="modelValue"
        :disabled="disabled"
        @input="$emit('update:modelValue', parseInt($event.target.value))"
      >
    </div>
  </div>
</template>

<script>
export default {
  props: {
    modelValue: {
      type: Number,
      default: 50
    },
    selectedGene: {
      type: String,
      required: true
    },
    label: {
      type: String,
      required: true
    },
    disabled: {
      type: Boolean,
      required: true
    }
  },
  emits: ["update:modelValue"]
};
</script>

<style>
.baseinput {
  width: var(--input-width);
  max-width: 90%;
  font-size: var(--fnt-normal);
  margin-bottom: 1em;
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.input {
  display: flex;
  flex-direction: column;
  align-items: center;
}

input {
  width: 8em;
}
</style>

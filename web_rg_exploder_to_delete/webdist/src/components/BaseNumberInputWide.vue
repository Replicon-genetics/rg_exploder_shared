<template>
  <div class="baseinput">
    <label v-if="label"> {{ label }} ({{ min }} to {{ max }}) </label>
    <input
      v-bind="$attrs"
      :value="modelValue"
      :min="min"
      :max="max"
      step="1"
      type="number"
      :disabled="disabled"
      @input="$emit('update:modelValue', Number($event.target.value))"
    />
  </div>
</template>

<script>
export default {
  props: {
    modelValue: {
      type: Number,
      default: 0
    },
    label: {
      type: String,
      default: ''
    },
    min: {
      type: String,
      default: '1'
    },
    max: {
      type: String,
      default: '5000'
    },
    disabled: {
      type: Boolean,
      required: true
    }
  },
  emits: ['update:modelValue'],
  watch: {
    modelValue: function() {
      console.log(`${this.label} is now ${this.modelValue}`)
      if (this.modelValue < this.min || this.modelValue > this.max) {
        const number = Math.min(this.max, Math.max(this.min, this.modelValue))
        this.$emit('update:modelValue', Number(number))
      }
    }
  }
}
</script>

<style scoped>
.baseinput {
  width: 100%;
  font-size: var(--fnt-normal);
  display: flex;
  flex-direction: row;
  justify-content: space-between;
}

input {
  text-align:right;
  width: 6em;
}
</style>

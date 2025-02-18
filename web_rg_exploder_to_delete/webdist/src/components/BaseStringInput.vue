<template>
  <div class="baseinput">
    <label v-if="label"> {{ label }} </label>
    <input
      v-bind="$attrs"
      :value="modelValue"
      type="string"
      :disabled="disabled"
      @input="$emit('update:modelValue', $event.target.value)"
    />
  </div>
</template>

<script>
export default {
  props: {
    modelValue: {
      type: String,
      default: 'AAAAAAA'
    },
    label: {
      type: String,
      default: ''
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
      // Some logic here to constrain string content?
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
  width: 20em;
  text-indent: 3px;
}
</style>

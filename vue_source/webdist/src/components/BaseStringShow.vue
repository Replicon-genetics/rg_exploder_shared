<template>
  <div class="baseinput">
    <label v-if="label"> {{ label }} </label>
    <input
      v-bind="$attrs"
      :value="modelValue"
      type="string"
      :disabled="disabled"
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
  //emits: ['update:modelValue'],
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
  width: 30em;
  text-indent: 3px;
}
</style>

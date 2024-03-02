<template>
  <div class="baseinput">
    <label v-if="label"></label>
    <input
      v-bind="$attrs"
      :value="modelValue"
      :disabled="disabled"
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
    disabled: {
      type: Boolean,
      required: true
    }
  },
  //emits: ['update:modelValue'],
  watch: {
    modelValue: function() {
      //console.log(`${this.label} is now ${this.modelValue}`)
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
  justify-content: left;
}

input {
  text-align:right;
  width: 7em;
}
</style>

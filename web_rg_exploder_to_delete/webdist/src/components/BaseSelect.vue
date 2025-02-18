<template>
  <div class="baseinput">
    <label>
      <a v-if="link" :href="link" target="_blank" class="link">{{ label }}</a>
      <template v-else>{{ label }}</template>
    </label>

    <select
      :disabled="disabled"
      :name="label"
      :value="modelValue"
      v-bind="{
        ...$attrs,
        onChange: $event => $emit('update:modelValue', $event.target.value)
      }"
    >
      <option
        v-for="option in options"
        :key="option"
        :value="option"
        :selected="option === modelValue"
      >
        {{ option }}
      </option>
    </select>
  </div>
</template>

<script>
export default {
  props: {
    modelValue: {
      type: [String, Number],
      default: ''
    },
    label: {
      type: String,
      default: ''
    },
    options: {
      type: Array,
      required: true
    },
    disabled: {
      type: Boolean,
      required: true
    },
    link: {
      type: String,
      default: ''
    }
  },
  emits: ['update:modelValue'],
  watch: {
    modelValue: function() {
      console.log(`${this.modelValue} is now selected`)
    }
  }
}
</script>

<style scoped>
.baseinput {
  max-width: 90%;
  font-size: var(--fnt-normal);
  margin-bottom: 1em;
  display: flex;
  justify-content: space-between;
}

select {
  width: 10em;
  text-indent: 3px;
}

.link {
  color: var(--primary-color);
}
</style>

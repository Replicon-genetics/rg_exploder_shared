<template>
  <div class="baseinput">
    <label v-if="label"> {{ label }} </label>
    <input
      v-bind="$attrs"
      :value="inputSequence"
      type="String"
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
      default: '-'
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
  data() {
    return {
      inputSequence: this.modelValue,
    };
  },
  methods: {
    handleInput() {
      // Restrict the input to accept only IUPAC DNA letters
      this.inputSequence = this.inputSequence.replace(/[^ACGTRYSWKMBDHVN]/gi, "");
      this.inputSequence = this.inputSequence.toUpperCase();
      // If empty, set to - (delete)
      if (this.inputSequence=="")
      {this.inputSequence="-"}

      // Emit the updated value
      this.$emit('update:modelValue', this.inputSequence);
    }
  },
  watch: {
    // Seemed to need this to set the inputSequence value to the incoming modelValue
    inputSequence: function() {
      this.handleInput();
    },
    modelValue: function(newValue) {
      this.inputSequence = newValue;
      //console.log(`${this.label} is now ${this.inputSequence}`)
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

<template>
<v-row no-gutters :style="{ padding: '0 10px' }">
    <v-col cols='4' :style="{ padding: '0 5px', margin: 0 }">
        <div class="baseinput">
            <label v-if="Blabel"> {{ Blabel }} {{"BARF" }} ({{ Bmin }} to {{ Bmax }}) </label>
            <input
            v-bind="$attrs"
            :value="BmodelValue"
            :Bmin="min"
            :Bmax="max"
            step="1"
            type="number"
            :disabled="disabled"
            @input="$emit('update:BmodelValue', Number($event.target.value))"
            />
        </div>
    </v-col>
    <v-col cols='5' :style="{ padding: '0 5px', margin: 0 }">
        <div class="baseinput">
            <label v-if="Elabel"> {{ Elabel }} {{"BARF2" }} ({{ Emin }} to {{ Emax }}) </label>
            <input
            v-bind="$attrs"
            :value="EmodelValue"
            :Emin="min"
            :Emax="max"
            step="1"
            type="number"
            :disabled="disabled"
            @input="$emit('update:EmodelValue', Number($event.target.value))"
            />
        </div>

    </v-col>
    <v-col cols='2' :style="{ padding: '0 5px', margin: 0 }">
        <v-btn>
            <v-icon>X</v-icon>
        </v-btn>
    </v-col>
</v-row>
  </template>
  
  <script>
  export default {
    props: {
      BmodelValue: {
        type: Number,
        default: 0
      },
      Blabel: {
        type: String,
        default: ''
      },
      Bmin: {
        type: String,
        default: '1'
      },
      Bmax: {
        type: String,
        default: '5000'
      },
      EmodelValue: {
        type: Number,
        default: 0
      },
      Elabel: {
        type: String,
        default: ''
      },
      Emin: {
        type: String,
        default: '1'
      },
      Emax: {
        type: String,
        default: '5000'
      },
      disabled: {
        type: Boolean,
        required: true
      }
    },
    emits: ['update:BmodelValue'],
    watch: {
      BmodelValue: function() {
        console.log(`${this.Blabel} is now ${this.BmodelValue}`)
        if (this.BmodelValue < this.Bmin || this.BmodelValue > this.Bmax) {
          const number = Math.min(this.Bmax, Math.max(this.Bmin, this.BmodelValue))
          this.$emit('update:BmodelValue', Number(number))
        }
      },
      EmodelValue: function() {
        console.log(`${this.Elabel} is now ${this.EmodelValue}`)
        if (this.EmodelValue < this.Emin || this.EmodelValue > this.Emax) {
          const number = Math.min(this.Emax, Math.max(this.Emin, this.EmodelValue))
          this.$emit('update:EmodelValue', Number(number))
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
  
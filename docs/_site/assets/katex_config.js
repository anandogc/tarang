{
    function renderKatex() {
        let macros = {}
        if (customElements) {
            class KatexInline extends HTMLElement {
                constructor() {
                    super();
                    katex.render(this.innerText, this, {throwOnError: false, displayMode: false, macros: macros, output: "html"});
                }
            }
            customElements.define("katex-inline", KatexInline)

            class KatexBlock extends HTMLElement {
                constructor() {
                    super();
                    katex.render(this.innerText, this, {throwOnError: true, displayMode: true, macros: macros, output: "html"});
                }
            }
            customElements.define("katex-block", KatexBlock)
        } else {
            document.querySelectorAll("katex-inline").forEach(
                (el) => {
                    katex.render(el.innerText, el, {throwOnError: false, displayMode: false, macros: macros, output: "html"});
                }
            )
            document.querySelectorAll("katex-block").forEach(
                (el) => {
                    katex.render(el.innerText, el, {throwOnError: false, displayMode: true, macros: macros, output: "html"});
                }
            )
        }
    }
}
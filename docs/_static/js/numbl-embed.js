/*
Define a custom HTML element <numbl-embed> that embeds a numbl iframe.

Usage:
<numbl-embed>
<iframe width="100%" height="500" frameborder="0"></iframe>

<script type="text/plain" class="numbl-script">
... numbl script ...
</script>

</numbl-embed>

Optional preamble: setup code that runs before the visible script on every Run
(e.g. installing a package with `mip load --install ...`). It is kept out of the
editor, and its console output is hidden behind a "Preparing…" message unless it
errors — in which case the output and error are shown. Override that message
per-embed with the `preparing-label` attribute (e.g. "Installing…").

<numbl-embed preparing-label="Installing…">
<iframe width="100%" height="500" frameborder="0"></iframe>

<script type="text/plain" class="numbl-preamble">
mip load --install owner/repo/package
</script>

<script type="text/plain" class="numbl-script">
... numbl script ...
</script>

</numbl-embed>

The classes `matlab-script` / `matlab-preamble` are accepted as synonyms for
`numbl-script` / `numbl-preamble` (backward compatibility with older embeds).

You can also specify the script URL via attribute (a `numbl-preamble` block still
applies):
<numbl-embed script="relative-or-absolute-url-to-script">
<iframe width="100%" height="500" frameborder="0"></iframe>
</numbl-embed>

Lazy mode (recommended when many embeds share a page, e.g. documentation):
add the `lazy` attribute and the iframe is not loaded until the reader clicks
an "Edit & run" button. This keeps the page light — no worker, editor, or
package download boots until the reader asks for it. Customize the button text
with the `label` attribute.

<numbl-embed lazy label="▶ Edit & run this example"
             script="https://example.com/snippet.m">
<iframe width="100%" height="560" frameborder="0"></iframe>
</numbl-embed>
*/
// Base64-encode a string as UTF-8. Plain btoa() throws on any character
// outside Latin-1 (e.g. an em dash or a Greek letter in a comment/title), so
// we go through the UTF-8 bytes. For pure-ASCII input this is byte-identical
// to btoa(), so it stays compatible with older /embed decoders.
function utf8ToBase64(text) {
  const bytes = new TextEncoder().encode(text);
  let binary = "";
  const chunk = 0x8000;
  for (let i = 0; i < bytes.length; i += chunk) {
    binary += String.fromCharCode.apply(null, bytes.subarray(i, i + chunk));
  }
  return btoa(binary);
}

class NumblEmbed extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", () => this.setup());
    } else {
      this.setup();
    }
  }

  setup() {
    this.iframe = this.querySelector("iframe");
    if (!this.iframe) {
      console.error("Missing iframe element in numbl-embed");
      return;
    }

    // Lazy mode: defer all loading until the reader opts in by clicking.
    if (this.hasAttribute("lazy")) {
      this.renderActivateButton();
    } else {
      this.activate();
    }
  }

  renderActivateButton() {
    this.iframe.style.display = "none";

    const button = document.createElement("button");
    button.type = "button";
    button.className = "numbl-embed-activate";
    button.textContent = this.getAttribute("label") || "▶ Edit & run in numbl";
    // Inline styles so the button looks reasonable even without external CSS;
    // a page may override via the .numbl-embed-activate class.
    button.style.cssText = [
      "display:inline-flex",
      "align-items:center",
      "gap:0.4em",
      "padding:0.45em 0.9em",
      "font-size:0.85rem",
      "font-family:inherit",
      "color:#fff",
      "background:#2e7d32",
      "border:none",
      "border-radius:4px",
      "cursor:pointer",
    ].join(";");

    button.addEventListener(
      "click",
      () => {
        button.remove();
        this.iframe.style.display = "";
        this.activate();
      },
      { once: true }
    );

    this.activateButton = button;
    this.insertBefore(button, this.iframe);
  }

  activate() {
    const defaultNumblUrl = "https://numbl.org";
    const numblUrl = this.attributes["numbl-url"]?.value || defaultNumblUrl;
    const cacheBust = `_cb=${Date.now()}`;
    const mode = this.getAttribute("mode");

    // REPL mode: no script needed, just load the REPL page
    if (mode === "repl") {
      this.iframe.src = `${numblUrl}/embed-repl?${cacheBust}`;
      return;
    }

    const encodePlain = text => {
      const text2 = text.replace(/&lt;/g, "<").replace(/&gt;/g, ">");
      return utf8ToBase64(text2);
    };

    // Optional preamble: setup code that runs (hidden) before the visible
    // script. `numbl-preamble` is preferred; `matlab-preamble` is a synonym.
    const preambleElement = this.querySelector(
      "script.numbl-preamble, script.matlab-preamble"
    );
    let preambleParam = preambleElement
      ? `&preamble=${encodePlain(preambleElement.textContent.trim())}`
      : "";

    // Message shown while the preamble runs (defaults to "Preparing..." in the
    // embed page). Override per-embed with the `preparing-label` attribute,
    // e.g. preparing-label="Installing...".
    const preparingLabel = this.getAttribute("preparing-label");
    if (preambleParam && preparingLabel) {
      preambleParam += `&preparing=${utf8ToBase64(preparingLabel)}`;
    }

    // `numbl-script` is preferred; `matlab-script` is a synonym.
    const scriptElement = this.querySelector(
      "script.numbl-script, script.matlab-script"
    );

    let scriptBase64;
    if (scriptElement) {
      scriptBase64 = encodePlain(scriptElement.textContent.trim());
    } else if (this.attributes.script) {
      const scriptUrl = this.attributes.script.value;
      this.loadScriptFromUrl(scriptUrl, preambleParam);
      return;
    } else {
      scriptBase64 = null;
    }

    if (scriptBase64) {
      this.iframe.src = `${numblUrl}/embed?script=${scriptBase64}${preambleParam}&${cacheBust}`;
    } else {
      this.iframe.src = `${numblUrl}/embed?${cacheBust}${preambleParam}`;
    }
  }

  async loadScriptFromUrl(url, preambleParam = "") {
    try {
      let absoluteUrl = url;
      if (url.startsWith("./") || url.startsWith("../")) {
        const baseUrl = window.location.href;
        absoluteUrl = new URL(url, baseUrl).href;
      }

      const response = await fetch(absoluteUrl);
      if (!response.ok) {
        throw new Error(`Failed to fetch script: ${response.statusText}`);
      }

      const scriptContent = await response.text();
      const scriptBase64 = utf8ToBase64(scriptContent);

      const defaultNumblUrl = "https://numbl.org";
      const numblUrl = this.attributes["numbl-url"]?.value || defaultNumblUrl;

      const cacheBust = `_cb=${Date.now()}`;
      this.iframe.src = `${numblUrl}/embed?script=${scriptBase64}${preambleParam}&${cacheBust}`;
    } catch (error) {
      console.error("Error loading script:", error);
      this.iframe.srcdoc = `
        <html>
          <body style="font-family: Arial; padding: 20px;">
            <h3>Error loading MATLAB script</h3>
            <p>${error.message}</p>
          </body>
        </html>
      `;
    }
  }
}

customElements.define("numbl-embed", NumblEmbed);

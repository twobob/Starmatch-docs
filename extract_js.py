"""
Extract JavaScript code from HTML documentation files.
This script parses Docco-style documentation HTML files and reconstructs the original JavaScript files.
"""

import re
import os
from html import unescape

# HTML files that document JavaScript files (from index.html)
html_to_js_mapping = {
    "coeffanalyse.html": "CoeffAnalyse.js",
    "engine.html": "engine.js",
    "exporter.html": "exporter.js",
    "global.html": "global.js",
    "importer.html": "importer.js",
    "index.html": "index.js",
    "locationpicker.jquery.html": "locationpicker.jquery.js",
    "output.html": "output.js",
    "privacy-policy.html": "privacy-policy.js",
    "record.html": "record.js",
    "support.html": "support.js"
}

def extract_code_from_html(html_file):
    """Extract JavaScript code from a Docco-style HTML documentation file."""
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find all code cells - they're in <td class="code"> elements
    # The code is within <pre><code class='prettyprint'>...</code></pre>
    pattern = r'<td class="code">\s*<pre><code class=\'prettyprint\'>(.*?)</code></pre>\s*</td>'
    matches = re.findall(pattern, content, re.DOTALL)
    
    # Combine all code snippets
    code_lines = []
    for match in matches:
        # Unescape HTML entities
        unescaped = unescape(match)
        # Only add non-empty code blocks
        if unescaped.strip():
            code_lines.append(unescaped)
    
    return ''.join(code_lines)

def main():
    """Extract JavaScript from all HTML documentation files."""
    script_dir = "script"
    
    # Ensure script directory exists
    os.makedirs(script_dir, exist_ok=True)
    
    # Process each HTML file
    for html_file, js_file in html_to_js_mapping.items():
        if os.path.exists(html_file):
            print(f"Processing {html_file} -> {js_file}")
            
            try:
                js_code = extract_code_from_html(html_file)
                
                # Write the JavaScript file
                output_path = os.path.join(script_dir, js_file)
                with open(output_path, 'w', encoding='utf-8') as f:
                    f.write(js_code)
                
                print(f"  ✓ Created {output_path} ({len(js_code)} characters)")
            except Exception as e:
                print(f"  ✗ Error processing {html_file}: {e}")
        else:
            print(f"  ⚠ {html_file} not found, skipping")
    
    print("\nExtraction complete!")

if __name__ == "__main__":
    main()

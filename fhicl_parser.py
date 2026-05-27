import re
import os

class FhiclParseException(Exception):
    pass

class FhiclParser:
    def __init__(self):
        self.definitions = {}  # The final configuration
        self.prolog = {}       # The PROLOG definitions
        self.in_prolog = False # State flag

    def parse_text(self, text):
        """
        Entry point: Parses a raw string.
        """
        clean_text = self._strip_comments(text)
        self._parse_document(clean_text, self.definitions)
        return self.definitions

    def _strip_comments(self, text):
        """
        Removes # and // comments.
        Ref: fhiclcpp/recursive_build_fhicl.hxx
        """
        lines = text.split('\n')
        cleaned = []
        for line in lines:
            # Remove // comments
            if "//" in line:
                line = line.split("//")[0]
            # Remove # comments
            if "#" in line:
                line = line.split("#")[0]
            cleaned.append(line)
        return "\n".join(cleaned)

    def _find_matching_bracket(self, text, start_idx, open_char, close_char):
        """
        Finds the index of the closing bracket balancing the opening one.
        Handles nested brackets and ignores brackets inside quotes.
        Ref: fhiclcpp/fhicl_doc.hxx :: find_matching_bracket
        """
        depth = 0
        in_quote = False
        quote_char = None
        
        for i in range(start_idx, len(text)):
            char = text[i]
            
            # Handle quotes to ignore brackets inside strings
            if char in ['"', "'"]:
                if not in_quote:
                    in_quote = True
                    quote_char = char
                elif char == quote_char:
                    in_quote = False
                continue
            
            if in_quote:
                continue

            if char == open_char:
                depth += 1
            elif char == close_char:
                depth -= 1
                if depth == 0:
                    return i
        
        raise FhiclParseException(f"Unmatched {open_char} starting near {text[start_idx:start_idx+10]}")

    def _parse_document(self, text, target_dict):
        """
        Main loop. Scans for Keys, Values, and Directives.
        Ref: fhiclcpp/recursive_build_fhicl.hxx :: parse_fhicl_document
        """
        idx = 0
        n = len(text)
        
        while idx < n:
            # Skip whitespace
            while idx < n and text[idx].isspace():
                idx += 1
            if idx >= n: break

            remaining = text[idx:]
            if remaining.startswith("BEGIN_PROLOG"):
                self.in_prolog = True
                idx += len("BEGIN_PROLOG")
                continue
            elif remaining.startswith("END_PROLOG"):
                self.in_prolog = False
                idx += len("END_PROLOG")
                continue

            # Identify Key
            colon_match = re.search(r':', text[idx:])
            if not colon_match:
                 break
            
            potential_key_end = idx + colon_match.start()
            
            # Check for directives (e.g. @table::)
            if text[idx] == '@':
                val, new_idx = self._parse_value(text, idx, target_dict if not self.in_prolog else self.prolog)
                if isinstance(val, dict):
                    (self.prolog if self.in_prolog else target_dict).update(val)
                idx = new_idx
                continue

            # Extract Key
            key_str = text[idx:potential_key_end].strip()
            idx = potential_key_end + 1 # Skip colon

            # Parse Value
            value, idx = self._parse_value(text, idx, target_dict if not self.in_prolog else self.prolog)

            # Insert into dictionary (handling dot notation a.b.c)
            current_level = self.prolog if self.in_prolog else target_dict
            
            parts = key_str.split('.')
            
            for i, part in enumerate(parts[:-1]):
                arr_match = re.match(r'(\w+)\[(\d+)\]', part)
                if arr_match:
                    pname, pidx = arr_match.groups()
                    pidx = int(pidx)
                    if pname not in current_level:
                        current_level[pname] = []
                    while len(current_level[pname]) <= pidx:
                        current_level[pname].append({})
                    current_level = current_level[pname][pidx]
                else:
                    if part not in current_level:
                        current_level[part] = {}
                    if not isinstance(current_level[part], dict):
                        current_level[part] = {}
                    current_level = current_level[part]

            final_key = parts[-1]
            arr_match = re.match(r'(\w+)\[(\d+)\]', final_key)
            if arr_match:
                pname, pidx = arr_match.groups()
                pidx = int(pidx)
                if pname not in current_level:
                    current_level[pname] = []
                while len(current_level[pname]) <= pidx:
                    current_level[pname].append(None) 
                current_level[pname][pidx] = value
            else:
                current_level[final_key] = value

    def _parse_value(self, text, idx, context):
        """
        Determines the type of value starting at idx and parses it.
        Ref: fhiclcpp/recursive_build_fhicl.hxx :: parse_object
        """
        # Skip whitespace
        while idx < len(text) and text[idx].isspace():
            idx += 1
        
        char = text[idx]

        # 1. Table: { ... }
        if char == '{':
            end_idx = self._find_matching_bracket(text, idx, '{', '}')
            content = text[idx+1 : end_idx]
            new_table = {}
            self._parse_document(content, new_table)
            return new_table, end_idx + 1

        # 2. Sequence: [ ... ]
        elif char == '[':
            end_idx = self._find_matching_bracket(text, idx, '[', ']')
            content = text[idx+1 : end_idx]
            elements = self._split_sequence(content)
            parsed_elements = []
            for el in elements:
                if not el.strip(): continue
                val, _ = self._parse_value(el, 0, context)
                parsed_elements.append(val)
            
            final_list = []
            for item in parsed_elements:
                if isinstance(item, dict) and "__splice_sequence__" in item:
                     final_list.extend(item["__splice_sequence__"])
                else:
                    final_list.append(item)
            
            return final_list, end_idx + 1

        # 3. Quoted Strings (Corrected)
        elif char == '"' or char == "'":
            # Simple search for the closing quote.
            # Ref: fhiclcpp/recursive_build_fhicl.hxx :: case '\'' and '\"'
            end_idx = text.find(char, idx + 1)
            
            if end_idx == -1:
                # C++ source also warns if strings span multiple lines, implying they shouldn't.
                raise FhiclParseException(f"Unmatched {char} starting near {text[idx:idx+20]}...")
            
            # Extract content without quotes
            val_str = text[idx+1 : end_idx]
            return val_str, end_idx + 1

        # 4. Directives: @...
        elif char == '@':
            match = re.search(r'[\s,}\]]', text[idx:])
            if match:
                end_token = idx + match.start()
            else:
                end_token = len(text)
            
            token = text[idx:end_token]
            
            if token == "@nil":
                return None, end_token

            if "::" in token:
                directive, target = token.split("::", 1)
                resolved_value = self._resolve_reference(target, context)
                
                if directive == "@local":
                    return resolved_value, end_token
                elif directive == "@table":
                    if not isinstance(resolved_value, dict):
                         raise FhiclParseException(f"@table directive target {target} is not a table")
                    return resolved_value, end_token
                elif directive == "@sequence":
                     if not isinstance(resolved_value, list):
                         raise FhiclParseException(f"@sequence directive target {target} is not a sequence")
                     return {"__splice_sequence__": resolved_value}, end_token
            
            return token, end_token

        # 5. Unquoted Atoms
        else:
            match = re.search(r'[\s,}\]]', text[idx:])
            if match:
                end_atom = idx + match.start()
            else:
                end_atom = len(text)
            
            raw_val = text[idx:end_atom]
            return self._convert_atom(raw_val), end_atom

    def _split_sequence(self, text):
        """
        Splits a comma-separated string, respecting brackets/quotes.
        Ref: fhiclcpp/fhicl_doc.hxx :: get_list_elements
        """
        elements = []
        depth_curly = 0
        depth_square = 0
        in_quote = False
        quote_char = None
        last_idx = 0
        
        for i, char in enumerate(text):
            if char in ['"', "'"]:
                if not in_quote:
                    in_quote = True
                    quote_char = char
                elif char == quote_char:
                    in_quote = False
            elif not in_quote:
                if char == '{': depth_curly += 1
                elif char == '}': depth_curly -= 1
                elif char == '[': depth_square += 1
                elif char == ']': depth_square -= 1
                elif char == ',' and depth_curly == 0 and depth_square == 0:
                    elements.append(text[last_idx:i].strip())
                    last_idx = i + 1
        
        elements.append(text[last_idx:].strip())
        return elements

    def _convert_atom(self, val_str):
        """
        Converts string to int/float/bool if possible.
        Ref: fhiclcpp/string_parsers/from_string.hxx
        """
        val_str = val_str.strip()
        if val_str == "true" or val_str == "True": return True
        if val_str == "false" or val_str == "False": return False
        
        try:
            return int(val_str)
        except ValueError:
            pass
        
        try:
            return float(val_str)
        except ValueError:
            pass
        
        return val_str

    def _resolve_reference(self, key, context):
        """
        Look up a key in Definitions, Prolog, or Context.
        Ref: fhiclcpp/recursive_build_fhicl.hxx :: deep_copy_resolved_reference_value
        """
        def lookup(d, k_str):
            parts = k_str.split('.')
            curr = d
            try:
                for p in parts:
                    if isinstance(curr, dict):
                        curr = curr[p]
                    else:
                        return None
                return curr
            except KeyError:
                return None

        val = lookup(self.definitions, key)
        if val is not None: return val
        
        val = lookup(self.prolog, key)
        if val is not None: return val
            
        val = lookup(context, key)
        if val is not None: return val

        raise FhiclParseException(f"Could not resolve reference: {key}")

if __name__ == "__main__":
    fcl_content = """
    BEGIN_PROLOG
    prolog_val: "I am from prolog"
    prolog_table: { p_key: p_val }
    END_PROLOG

    group: {
      val1: 5
      val2: "string"
      nested: {
        a: 1
        b: [1, 2, 3]
      }
    }
    
    my_ref: @local::group.nested.a
    from_prolog: @local::prolog_val
    """
    

    parser = FhiclParser()
    try:
        result = parser.parse_text(fcl_content)
        import json
        print(json.dumps(result, indent=2))
    except FhiclParseException as e:
        print(f"Error: {e}")
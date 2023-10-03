#include "xll_FRE6233.h"

using namespace xll;

// Create XML documentation and index.html in `$(TargetPath)` folder.
// Use `xsltproc file.xml -o file.html` to create HTML documentation.
#ifdef _DEBUG

xll_url_set FRE6233("https://keithalewis.github.io/FRE6233/");
Auto<Open> xao_template_docs([]() {

	return Documentation(CATEGORY, "Documentation for " CATEGORY ".");

});

#endif // _DEBUG

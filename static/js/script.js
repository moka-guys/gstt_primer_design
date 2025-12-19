let igvBrowser = null;

document.addEventListener("DOMContentLoaded", async function () {
    const igvDiv = document.getElementById("igv-div");
    if (!igvDiv) return;

    const options = {
        genome: initialQuery.genome,
        showReference: false,
        locus: initialQuery.locus,
        tracks: [
            {
                name: "Primers",
                type: "annotation",
                format: "bed",
                url: primersBedUrl,
                displayMode: "EXPANDED",
                color: "green"
            },
            {
                name: "SNPs",
                type: "annotation",
                format: "bed",
                url: snpsbedUrl,
                displayMode: "EXPANDED",
                color: "red"
            }
        ]
    };

    igvBrowser = await igv.createBrowser(igvDiv, options);
    console.log("IGV Reloaded with fresh BED files");
});

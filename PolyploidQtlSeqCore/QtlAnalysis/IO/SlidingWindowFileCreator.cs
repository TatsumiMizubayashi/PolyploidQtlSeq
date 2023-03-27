using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.IO
{
    /// <summary>
    /// Sliding Windowファイルクリエーター
    /// </summary>
    internal static class SlidingWindowFileCreator
    {
        internal const string QTL = SnpIndexFileCreator.QTL;
        private const string FILE_NAME = "SlidingWindow.txt";

        private static readonly string[] _fieldNames = new[]
            {
                "Chr",
                "Start",
                "End",
                "Center",
                "P99",
                "P95",
                "Variant Count",
                "P99 QTL Variant Count",
                "P99 +Delta SNP-index QTL Variant Count",
                "P99 -Delta SNP-index QTL Variant Count",
                "P99 QTL Variant Rate",
                "P95 QTL Variant Count",
                "P95 +Delta SNP-index QTL Variant Count",
                "P95 -Delta SNP-index QTL Variant Count",
                "P95 QTL Variant Rate",
                "Mean P99 Threshold",
                "Mean P95 Threshold",
                "Mean Bulk1 SNP-index",
                "Mean Bulk2 SNP-index",
                "Mean Delta SNP-index",
                "Mean P-Value",
                "Mean -log10(P)",
                "Bulk1 SNP-index=0 Variant Count",
                "Bulk2 SNP-index=0 Variant Count"
            };

        /// <summary>
        /// SlidingWindowファイルを作成する。
        /// </summary>
        /// <param name="outDir">ファイル保存ディレクトリ</param>
        /// <param name="windows">Window</param>
        public static void Create(OutputDirectory outDir, Window[] windows)
        {
            var filePath = outDir.CreateFilePath(FILE_NAME);
            using var writer = new StreamWriter(filePath);
            var fieldLine = string.Join("\t", _fieldNames);
            writer.WriteLine(fieldLine);

            foreach (var window in windows)
            {
                var values = window.VariantCount.Count == 0
                    ? GetEmptyWindowValues(window)
                    : GetWindowValues(window);
                var line = string.Join("\t", values);
                writer.WriteLine(line);
            }
        }


        /// <summary>
        /// 変異0個Windowの値を取得する。
        /// </summary>
        /// <param name="window">window</param>
        /// <returns>値</returns>
        private static string[] GetEmptyWindowValues(Window window)
        {
            return new[]
            {
                window.GenomePosition.ChrName,
                window.GenomePosition.Start.ToString(),
                window.GenomePosition.End.ToString(),
                window.GenomePosition.Center.ToString(),
                "",     // P99
                "",     // P95
                "0",    // Variant Count
                "",     // P99 QTL Variant Count
                "",     // P99 +ΔSNP-index QTL Variant Count
                "",     // P99 -ΔSNP-index QTL Variant Count
                "",     // P99 QTL Variant Rate
                "",     // P95 QTL Variant Count
                "",     // P95 +ΔSNP-index QTL Variant Count
                "",     // P95 -ΔSNP-index QTL Variant Count
                "",     // P95 QTL Variant Rate
                "",     // Mean P99 Threshold
                "",     // Mean P95 Threshold
                "",     // Mean Bulk1 SNP-index
                "",     // Mean Bulk2 SNP-index
                "",     // Mean ΔSNP-index
                "",     // Mean P-Value
                "",     // Mean -log10(P)
                "",     // Bulk1 SNP-index=0 Variant Count
                "",     // Bulk2 SNP-index=0 Variant Count
            };
        }

        /// <summary>
        /// Windowの値を取得する。
        /// </summary>
        /// <param name="window">window</param>
        /// <returns>値</returns>
        private static string[] GetWindowValues(Window window)
        {
            var genomePos = window.GenomePosition;
            var p99 = window.P99Qtl;
            var p95 = window.P95Qtl;
            var variantCount = window.VariantCount;

            return new[]
            {
                genomePos.ChrName,
                genomePos.Start.ToString(),
                genomePos.End.ToString(),
                genomePos.Center.ToString(),
                p99.IsQtl ? QTL : "",
                p95.IsQtl ? QTL : "",
                variantCount.Count.ToString(),
                p99.TotalQtlVariantCount.ToString(),
                p99.PlusDeltaSnpIndexQtlVariantCount.ToString(),
                p99.MinusDeltaSnpIndexQtlVariantCount.ToString(),
                p99.TotalQtlVariantRate.ToString(),
                p95.TotalQtlVariantCount.ToString(),
                p95.PlusDeltaSnpIndexQtlVariantCount.ToString(),
                p95.MinusDeltaSnpIndexQtlVariantCount.ToString(),
                p95.TotalQtlVariantRate.ToString(),
                p99.ThresholdDeltaSnpIndex.Value.ToString(),
                p95.ThresholdDeltaSnpIndex.Value.ToString(),
                window.AverageBulk1SnpIndex.Value.ToString(),
                window.AverageBulk2SnpIndex.Value.ToString(),
                window.AverageDeltaSnpIndex.Value.ToString(),
                window.AveragePValue.Value.ToString(),
                window.AverageScore.Value.ToString(),
                variantCount.Bulk1ZeroSnpIndexVariantCount.ToString(),
                variantCount.Bulk2ZeroSnpIndexVariantCount.ToString()
            };
        }

    }
}

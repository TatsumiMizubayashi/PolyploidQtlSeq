using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// Window
    /// </summary>
    internal class Window : IHasGenomePositionItem
    {
        /// <summary>
        /// Window情報を作成する。
        /// </summary>
        /// <param name="genomePosition">Window位置情報</param>
        /// <param name="averageBulk1SnpIndex">Bulk1 SNP-index平均値</param>
        /// <param name="averageBulk2SnpIndex">Bulk2 SNP-index平均値</param>
        /// <param name="averageDeltaSnpIndex">ΔSNP-index平均値</param>
        /// <param name="averageScore">スコア平均値</param>
        /// <param name="p95Qtl">P95QTL情報</param>
        /// <param name="p99Qtl">P99QTL情報</param>
        /// <param name="variantCount">変異数統計情報</param>
        public Window(GenomePosition genomePosition, Bulk1SnpIndex averageBulk1SnpIndex, Bulk2SnpIndex averageBulk2SnpIndex,
            DeltaSnpIndex averageDeltaSnpIndex, Score averageScore, WindowQtl p95Qtl, WindowQtl p99Qtl,
            VariantCount variantCount)
        {
            GenomePosition = genomePosition;
            AverageBulk1SnpIndex = averageBulk1SnpIndex;
            AverageBulk2SnpIndex = averageBulk2SnpIndex;
            AverageDeltaSnpIndex = averageDeltaSnpIndex;
            AverageScore = averageScore;
            AveragePValue = averageScore.ToPValue();
            P95Qtl = p95Qtl;
            P99Qtl = p99Qtl;
            VariantCount = variantCount;
        }


        /// <summary>
        /// Window位置を取得する。
        /// </summary>
        public GenomePosition GenomePosition { get; }

        /// <summary>
        /// Bulk1のSNP-index平均値を取得する。
        /// </summary>
        public Bulk1SnpIndex AverageBulk1SnpIndex { get; }

        /// <summary>
        /// Bulk2のSNP-index平均値を取得する。
        /// </summary>
        public Bulk2SnpIndex AverageBulk2SnpIndex { get; }

        /// <summary>
        /// ΔSNP-index平均値を取得する。
        /// </summary>
        public DeltaSnpIndex AverageDeltaSnpIndex { get; }

        /// <summary>
        /// スコア平均値を取得する。
        /// </summary>
        public Score AverageScore { get; }

        /// <summary>
        /// PValue平均値を取得する。
        /// (スコアから算出した値）
        /// </summary>
        public PValue AveragePValue { get; }

        /// <summary>
        /// P95QTL情報を取得する。
        /// </summary>
        public WindowQtl P95Qtl { get; }

        /// <summary>
        /// P99QTL情報を取得する。
        /// </summary>
        public WindowQtl P99Qtl { get; }

        /// <summary>
        /// 変異数統計情報を取得する。
        /// </summary>
        public VariantCount VariantCount { get; }

    }
}

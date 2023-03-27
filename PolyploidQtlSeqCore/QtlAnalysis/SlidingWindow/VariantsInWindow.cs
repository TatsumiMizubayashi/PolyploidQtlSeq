using Sequence.Position;

namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// Window内の変異
    /// </summary>
    internal class VariantsInWindow
    {
        private readonly GenomePosition _windowPosition;
        private readonly SnpIndexVariantWithQtl[] _variants;

        /// <summary>
        /// Window内の変異情報を作成する。
        /// </summary>
        /// <param name="windowPosition">window位置</param>
        /// <param name="variants">window内の変異</param>
        public VariantsInWindow(GenomePosition windowPosition, SnpIndexVariantWithQtl[] variants)
        {
            _windowPosition = windowPosition;
            _variants = variants;
        }

        /// <summary>
        /// Window情報に変換する。
        /// </summary>
        /// <returns>Window情報</returns>
        public Window ToWindow()
        {
            if (_variants.Length == 0) return CreateEmptyWindow();

            var averageBulk1SnpIndex = CalcAverageBulk1SnpIndex();
            var averageBulk2SnpIndex = CalcAverageBulk2SnpIndex();
            var averageDeltaSnpIndex = CalcAverageDeltaSnpIndex();
            var averageScore = CalcAverageScore();
            var p95WindowQtl = CreateP95WindowQtl(averageDeltaSnpIndex);
            var p99WindowQtl = CreateP99WindowQtl(averageDeltaSnpIndex);
            var variantCount = CreateVariantCount();

            return new Window(
                _windowPosition,
                averageBulk1SnpIndex,
                averageBulk2SnpIndex,
                averageDeltaSnpIndex,
                averageScore,
                p95WindowQtl,
                p99WindowQtl,
                variantCount);
        }

        /// <summary>
        /// Bulk1 SNP-index平均を計算する。
        /// </summary>
        /// <returns>Bulk1 SNP-index平均</returns>
        private Bulk1SnpIndex CalcAverageBulk1SnpIndex()
        {
            var average = _variants.Average(x => x.Bulk1.SnpIndex.Value);
            return new Bulk1SnpIndex(average);
        }

        /// <summary>
        /// Bulk2 SNP-index平均を計算する。
        /// </summary>
        /// <returns>Bulk2 SNP-index平均</returns>
        private Bulk2SnpIndex CalcAverageBulk2SnpIndex()
        {
            var average = _variants.Average(x => x.Bulk2.SnpIndex.Value);
            return new Bulk2SnpIndex(average);
        }

        /// <summary>
        /// ΔSNP-index平均を計算する。
        /// </summary>
        /// <returns>ΔSNP-index平均</returns>
        private DeltaSnpIndex CalcAverageDeltaSnpIndex()
        {
            var average = _variants.Average(x => x.DeltaSnpIndex.Value);
            return new DeltaSnpIndex(average);
        }

        /// <summary>
        /// スコア平均を計算する。
        /// </summary>
        /// <returns>スコア平均</returns>
        private Score CalcAverageScore()
        {
            var average = _variants.Average(x => x.Score.Value);
            return new Score(average);
        }

        /// <summary>
        /// P95ΔSNP-indexしきい値平均を計算する。
        /// </summary>
        /// <returns>P95ΔSNP-indexしきい値平均</returns>
        private QtlThresholdDeltaSnpIndex CalcAverageP95ThresholdDeltaSnpIndex()
        {
            var average = _variants.Average(x => x.P95ThresholdDeltaSnpIndex.Value);
            return new QtlThresholdDeltaSnpIndex(average);
        }

        /// <summary>
        /// P99ΔSNP-indexしきい値平均を計算する。
        /// </summary>
        /// <returns>P99ΔSNP-indexしきい値平均</returns>
        private QtlThresholdDeltaSnpIndex CalcAverageP99ThresholdDeltaSnpIndex()
        {
            var average = _variants.Average(x => x.P99ThresholdDeltaSnpIndex.Value);
            return new QtlThresholdDeltaSnpIndex(average);
        }

        /// <summary>
        /// P95WindowQTL情報を作成する。
        /// </summary>
        /// <param name="deltaSnpIndex">ΔSNP-index</param>
        /// <returns>P95WindowQTL情報</returns>
        private WindowQtl CreateP95WindowQtl(DeltaSnpIndex deltaSnpIndex)
        {
            var threshold = CalcAverageP95ThresholdDeltaSnpIndex();
            var qtlVariants = _variants.Where(x => x.P95Qtl).ToArray();
            var plusQtlCount = GetPlusDeltaSnpIndexCount(qtlVariants);
            var minusQtlCount = GetMinusDeltaSnpIndexCount(qtlVariants);
            var isQtl = deltaSnpIndex.IsQtl(threshold);

            return new WindowQtl(threshold, isQtl, plusQtlCount, minusQtlCount, _variants.Length);
        }

        /// <summary>
        /// P99WindowQTL情報を作成する。
        /// </summary>
        /// <param name="deltaSnpIndex">ΔSNP-index</param>
        /// <returns>P99WindowQTL情報</returns>
        private WindowQtl CreateP99WindowQtl(DeltaSnpIndex deltaSnpIndex)
        {
            var threshold = CalcAverageP99ThresholdDeltaSnpIndex();
            var qtlVariants = _variants.Where(x => x.P99Qtl).ToArray();
            var plusQtlCount = GetPlusDeltaSnpIndexCount(qtlVariants);
            var minusQtlCount = GetMinusDeltaSnpIndexCount(qtlVariants);
            var isQtl = deltaSnpIndex.IsQtl(threshold);

            return new WindowQtl(threshold, isQtl, plusQtlCount, minusQtlCount, _variants.Length);
        }

        /// <summary>
        /// 変異数統計情報を作成する。
        /// </summary>
        /// <returns>変異数統計情報</returns>
        private VariantCount CreateVariantCount()
        {
            var bulk1ZeroSnpIndexCount = _variants.Count(x => x.Bulk1.SnpIndex.Value == 0);
            var bulk2ZeroSnpIndexCount = _variants.Count(x => x.Bulk2.SnpIndex.Value == 0);

            return new VariantCount(_variants.Length, bulk1ZeroSnpIndexCount, bulk2ZeroSnpIndexCount);
        }

        /// <summary>
        /// ΔSNP-indexが正の値を持つ変異数を取得する。
        /// </summary>
        /// <param name="qtlVariants">QTL変異</param>
        /// <returns>ΔSNP-indexが正の値を持つ変異数</returns>
        private static int GetPlusDeltaSnpIndexCount(SnpIndexVariantWithQtl[] qtlVariants)
        {
            return qtlVariants.Count(x => x.DeltaSnpIndex.Value > 0);
        }

        /// <summary>
        /// ΔSNP-indexが負の値を持つ変異数を取得する。
        /// </summary>
        /// <param name="qtlVariants">QTL変異</param>
        /// <returns>ΔSNP-indexが負の値を持つ変異数</returns>
        private static int GetMinusDeltaSnpIndexCount(SnpIndexVariantWithQtl[] qtlVariants)
        {
            return qtlVariants.Count(x => x.DeltaSnpIndex.Value < 0);
        }


        /// <summary>
        /// 変異0のWindowを作成する。
        /// </summary>
        /// <returns>Window</returns>
        private Window CreateEmptyWindow()
        {
            return new Window(
                _windowPosition,
                new Bulk1SnpIndex(0),
                new Bulk2SnpIndex(0),
                new DeltaSnpIndex(0),
                new Score(0),
                new WindowQtl(new QtlThresholdDeltaSnpIndex(0), isQtl: false, 0, 0, 0),
                new WindowQtl(new QtlThresholdDeltaSnpIndex(0), isQtl: false, 0, 0, 0),
                CreateVariantCount());
        }
    }
}

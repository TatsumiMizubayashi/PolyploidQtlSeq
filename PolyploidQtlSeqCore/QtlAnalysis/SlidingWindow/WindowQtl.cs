namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// WindowのQTL情報
    /// </summary>
    internal class WindowQtl
    {
        /// <summary>
        /// WindowのQTL情報を作成する。
        /// </summary>
        /// <param name="threshold">QTLとみなすΔSNPindexのしきい値</param>
        /// <param name="isQtl">QTLかどうか</param>
        /// <param name="plusDeltaSnpIndexQtlCount">ΔSNP-indexが正の値のQTL変異数</param>
        /// <param name="minusDeltaSnpIndexQtlCount">ΔSNP-indexが負の値のQTL変異数</param>
        /// <param name="variantCount">window内の変異数</param>
        public WindowQtl(QtlThresholdDeltaSnpIndex threshold, bool isQtl, int plusDeltaSnpIndexQtlCount, int minusDeltaSnpIndexQtlCount, int variantCount)
        {
            if (plusDeltaSnpIndexQtlCount < 0) throw new ArgumentException(null, nameof(plusDeltaSnpIndexQtlCount));
            if (minusDeltaSnpIndexQtlCount < 0) throw new ArgumentException(null, nameof(minusDeltaSnpIndexQtlCount));
            if (variantCount < 0) throw new ArgumentException(null, nameof(variantCount));

            ThresholdDeltaSnpIndex = threshold;
            IsQtl = isQtl;
            PlusDeltaSnpIndexQtlVariantCount = plusDeltaSnpIndexQtlCount;
            MinusDeltaSnpIndexQtlVariantCount = minusDeltaSnpIndexQtlCount;
            TotalQtlVariantCount = plusDeltaSnpIndexQtlCount + minusDeltaSnpIndexQtlCount;

            if (variantCount == 0)
            {
                PlusDeltaSnpIndexQtlVariantRate = 0;
                MinusDeltaSnpIndexQtlVariantRate = 0;
                TotalQtlVariantRate = 0;
            }
            else
            {
                var total = (double)variantCount;
                PlusDeltaSnpIndexQtlVariantRate = PlusDeltaSnpIndexQtlVariantCount / total;
                MinusDeltaSnpIndexQtlVariantRate = MinusDeltaSnpIndexQtlVariantCount / total;
                TotalQtlVariantRate = TotalQtlVariantCount / total;
            }
        }

        /// <summary>
        /// QTLとみなすΔSNP-indexしきい値を取得する。
        /// </summary>
        public QtlThresholdDeltaSnpIndex ThresholdDeltaSnpIndex { get; }

        /// <summary>
        /// ΔSNP-indexが正の数のQTL変異数を取得する。
        /// </summary>
        public int PlusDeltaSnpIndexQtlVariantCount { get; }

        /// <summary>
        /// ΔSNP-indexが正の数のQTL変異割合を取得する。
        /// </summary>
        public double PlusDeltaSnpIndexQtlVariantRate { get; }

        /// <summary>
        /// ΔSNP-indexが負の数のQTL変異数を取得する。
        /// </summary>
        public int MinusDeltaSnpIndexQtlVariantCount { get; }

        /// <summary>
        /// ΔSNP-indexが負の数のQTL変異割合を取得する。
        /// </summary>
        public double MinusDeltaSnpIndexQtlVariantRate { get; }

        /// <summary>
        /// QTL変異数を取得する。
        /// </summary>
        public int TotalQtlVariantCount { get; }

        /// <summary>
        /// QTL変異割合を取得する。
        /// </summary>
        public double TotalQtlVariantRate { get; }

        /// <summary>
        /// QTLかどうかを取得する。
        /// </summary>
        public bool IsQtl { get; }
    }
}

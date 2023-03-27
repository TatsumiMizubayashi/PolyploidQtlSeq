namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// F1 QTLなしΔSNP-index分布
    /// </summary>
    internal class F1NoQtlDeltaSnpIndexDistribution
    {
        /// <summary>
        /// 指定値以上になるΔSNP-index値が0個の場合の代替値
        /// </summary>
        private const double ALT_COUNT_ZERO = 0.1;

        private readonly double[] _sortedAbsDeltaSnpIndice;

        /// <summary>
        /// F1 QTLなしΔSNP-index分布を作成する。
        /// </summary>
        /// <param name="deltaSnpIndice">ΔSNP-indice</param>
        public F1NoQtlDeltaSnpIndexDistribution(IEnumerable<DeltaSnpIndex> deltaSnpIndice)
        {
            _sortedAbsDeltaSnpIndice = deltaSnpIndice
                .Select(x => x.Abs())
                .OrderBy(x => x)
                .ToArray();

            if (_sortedAbsDeltaSnpIndice.Length == 0) throw new ArgumentException("ΔSNP-index is 0.");
        }

        /// <summary>
        /// 指定した優位水準に該当するしきい値を取得する。
        /// </summary>
        /// <param name="significanceLevel">優位水準</param>
        /// <returns>しきい値</returns>
        public QtlThresholdDeltaSnpIndex GetThreshold(double significanceLevel)
        {
            if (significanceLevel < 0 || significanceLevel > 1) throw new ArgumentException(null, nameof(significanceLevel));

            var index = (int)(_sortedAbsDeltaSnpIndice.Length * significanceLevel) - 1;

            return new QtlThresholdDeltaSnpIndex(_sortedAbsDeltaSnpIndice[index]);
        }

        /// <summary>
        /// 指定したΔSNP-indexのPValueを取得する。
        /// </summary>
        /// <param name="deltaSnpIndex">ΔSNP-index</param>
        /// <returns>P value</returns>
        public PValue GetPValue(DeltaSnpIndex deltaSnpIndex)
        {
            var absDeltaSnpIndex = deltaSnpIndex.Abs();

            var count = _sortedAbsDeltaSnpIndice.Count(x => x >= absDeltaSnpIndex);

            // P値が0だとスコア計算ができないので個数が0の場合は変わりの個数を使用してP値を計算する
            var pValue = count == 0
                ? ALT_COUNT_ZERO / _sortedAbsDeltaSnpIndice.Length
                : (double)count / _sortedAbsDeltaSnpIndice.Length;

            return new PValue(pValue);
        }
    }
}

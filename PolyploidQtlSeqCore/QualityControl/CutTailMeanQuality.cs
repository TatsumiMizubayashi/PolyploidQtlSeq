namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 3'末端トリム時の平均クオリティ
    /// </summary>
    internal class CutTailMeanQuality
    {
        /// <summary>
        /// 3'末端トリム時の平均クオリティの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 3'末端トリム時の平均クオリティの最大値
        /// </summary>
        private const int MAXIMUM = 30;


        /// <summary>
        /// 3'末端トリム時の平均クオリティを作成する。
        /// </summary>
        /// <param name="quality">平均クオリティ</param>
        public CutTailMeanQuality(int quality)
        {
            if (quality < MINIMUM || quality > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(quality));

            Value = quality;
        }

        /// <summary>
        /// 3'末端トリム時の平均クオリティを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"--cut_tail_mean_quality {Value}";
        }
    }
}

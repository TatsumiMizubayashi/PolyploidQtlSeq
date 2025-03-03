namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 塩基クオリティ
    /// </summary>
    internal class BaseQuality
    {
        /// <summary>
        /// 塩基クオリティの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 塩基クオリティの最大値
        /// </summary>
        private const int MAXIMUM = 50;

        /// <summary>
        /// 塩基クオリティを作成する。
        /// </summary>
        /// <param name="quality">塩基クオリティ</param>
        public BaseQuality(int quality)
        {
            if (quality < MINIMUM || quality > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(quality));

            Value = quality;
        }

        /// <summary>
        /// 塩基クオリティを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"-q {Value}";
        }
    }
}

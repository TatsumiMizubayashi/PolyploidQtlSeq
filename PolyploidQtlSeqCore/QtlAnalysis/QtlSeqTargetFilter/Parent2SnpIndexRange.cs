namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Parent2 SNP-index範囲
    /// </summary>
    internal class Parent2SnpIndexRange
    {
        private static readonly char[] _splitter = ['-'];

        /// <summary>
        /// SNP-indexの最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// SNP-indexの最大値
        /// </summary>
        private const double MAXIMUM = 1.0;


        /// <summary>
        /// Parent2 SNP-index範囲
        /// </summary>
        /// <param name="range">SNP-index範囲</param>
        public Parent2SnpIndexRange(string range)
        {
            if (string.IsNullOrWhiteSpace(range)) throw new ArgumentException(null, nameof(range));

            var items = range.Split(_splitter);
            if (items.Length != 2) throw new ArgumentException("Specify by lower limit - upper limit.", nameof(range));

            if (!double.TryParse(items[0], out var lower))
                throw new ArgumentException("The lower limit connot bo converted to a numerical value.", nameof(range));
            if (!double.TryParse(items[1], out var upper))
                throw new ArgumentException("The upper limit cannot be converted to a numerical value.", nameof(range));

            if (lower < MINIMUM || lower > MAXIMUM)
                throw new ArgumentException("The lower limit should be specified in the range of 0.0 to 1.0.", nameof(range));
            if (upper < MINIMUM || upper > MAXIMUM)
                throw new ArgumentException("The upper limit should be specified in the range of 0.0 to 1.0.", nameof(range));
            if (lower == upper) throw new ArgumentException("The lower and upper limits are the same value.", nameof(range));
            if (lower > upper) throw new ArgumentException("The lower limit is greater than the upper limit.", nameof(range));

            Lower = lower;
            Upper = upper;
        }

        /// <summary>
        /// 範囲の下限値を取得する。
        /// </summary>
        internal double Lower { get; }

        /// <summary>
        /// 範囲の上限値を取得する。
        /// </summary>
        internal double Upper { get; }
    }
}

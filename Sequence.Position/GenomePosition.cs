
namespace Sequence.Position
{
    /// <summary>
    /// ゲノム上の位置
    /// </summary>
    public sealed class GenomePosition
    {
        /// <summary>
        /// 位置情報が空のゲノム位置を作成する。
        /// </summary>
        public GenomePosition()
        {
            ChrName = "";
            Start = 0;
            End = 0;
            Length = 0;
            Center = 0;
            IsEmpty = true;
        }

        /// <summary>
        /// ゲノム位置を作成する。
        /// </summary>
        /// <param name="chr">染色体名</param>
        /// <param name="start">開始位置</param>
        /// <param name="end">終了位置</param>
        public GenomePosition(string chr, int start, int end)
        {
            if (string.IsNullOrWhiteSpace(chr)) throw new ArgumentException(null, nameof(chr));
            if (start < 1) throw new ArgumentException(null, nameof(start));
            if (start > end) throw new ArgumentException(null, nameof(end));

            ChrName = chr;
            Start = start;
            End = end;
            Length = end - start + 1;
            Center = start + Length / 2;
            IsEmpty = false;
        }

        /// <summary>
        /// 染色体名を取得する。
        /// </summary>
        public string ChrName { get; }

        /// <summary>
        /// 開始位置を取得する。
        /// </summary>
        public int Start { get; }

        /// <summary>
        /// 終了位置を取得する。
        /// </summary>
        public int End { get; }

        /// <summary>
        /// 領域長を取得する。
        /// </summary>
        public int Length { get; }

        /// <summary>
        /// 中心位置を取得する。
        /// </summary>
        public int Center { get; }

        /// <summary>
        /// 位置情報が空かどうかを取得する。
        /// </summary>
        public bool IsEmpty { get; }

        /// <summary>
        /// 位置が完全一致するか調査する。
        /// </summary>
        /// <param name="position">調査対象</param>
        /// <returns>完全一致するならtrue</returns>
        public bool IsMatch(GenomePosition position)
        {
            if (IsEmpty || position.IsEmpty) return false;

            if (ChrName != position.ChrName) return false;
            if (Start != position.Start) return false;
            if (End != position.End) return false;

            return true;
        }

        /// <summary>
        /// 位置に重なりがあるかを調査する。
        /// </summary>
        /// <param name="position">調査対象</param>
        /// <returns>重なりがあるならtrue</returns>
        public bool IsOverlap(GenomePosition position)
        {
            if (IsEmpty || position.IsEmpty) return false;

            if (ChrName != position.ChrName) return false;

            // 起点となる位置情報内に他位置のStartまたはEndがあれば重なっている
            return Start <= position.Start && position.Start <= End
                || Start <= position.End && position.End <= End
                || position.Start <= Start && Start <= position.End
                || position.Start <= End && End <= position.End;
        }


        /// <summary>
        /// 連結した位置情報を作成する。
        /// 2つの位置に重なりがある場合のみ連結できる。
        /// </summary>
        /// <param name="position">連結対象</param>
        /// <returns>連結後の位置情報</returns>
        public GenomePosition Merge(GenomePosition position)
        {
            if (IsEmpty || position.IsEmpty) throw new InvalidOperationException("位置情報が空のインスタンスはMergeできません。");
            if (!IsOverlap(position)) throw new ArgumentException("重なりがない位置情報はMergeできません。");

            var start = Math.Min(Start, position.Start);
            var end = Math.Max(End, position.End);

            return new GenomePosition(ChrName, start, end);
        }

        /// <summary>
        /// 拡張した位置情報を作成する。
        /// ±lengthした位置が作成される。
        /// </summary>
        /// <param name="length">拡張する長さ</param>
        /// <returns>拡張後の位置情報</returns>
        public GenomePosition Enlagement(int length)
        {
            if (IsEmpty) throw new InvalidOperationException("位置情報が空のインスタンスは拡張できません。");
            if (length < 1) throw new ArgumentException(null, nameof(length));

            var start = Start - length;
            if (start < 1) start = 1;

            var end = End + length;

            return new GenomePosition(ChrName, start, end);
        }
    }
}
